"""knobhole_lib.py -- build, simulate and score knob-peptide variants bound to
fibrinogen gamma-module hole 'a' (OpenMM, implicit solvent, single-trajectory
MM/GBSA).  CPU-only so the study can be reproduced without a GPU.

Model
-----
* Template: PDB 2FIB (human fibrinogen gamma-module, residues 143-411, 2.01 A,
  bound to GPRP; native gamma numbering).  Crystal waters are discarded.
* Pocket model: every receptor residue with a heavy atom within ``shell``
  angstrom of the template GPRP, plus the structural Ca2+ and its ligands when
  any of them falls in the shell.  Segment ends are capped with ACE/NME built
  from the real neighbouring atoms, so cuts carry no spurious charges.
  The shell is chosen once from the template peptide, so it is identical for
  every variant being compared.
* All receptor C-alpha atoms are restrained (harmonic, default 1 kcal/mol/A^2),
  side chains and the peptide are free.  A flat-bottom "leash" keeps the peptide
  backbone centroid within ``leash`` A of its starting position.
* Peptides are N-terminally free (NH3+, the new N-terminus after thrombin
  cleavage) and C-terminally capped with N-methylamide (NME) to mimic a
  continuing peptide bond rather than a free carboxylate.
* Positions are written P1..P4 (P1 = Gly of native GPRP).
* Solvent: OBC-II generalized Born (native OpenMM GBSAOBCForce, mbondi2 radii,
  solute/solvent dielectric 1/78.5, ACE non-polar term); no explicit salt term.
* Force field: AMBER ff14SB (+ TIP3P-set Ca2+).
"""
from __future__ import annotations

import math
import os
import time

import numpy as np
import openmm as mm
from openmm import app, unit
from openmm.app.internal.customgbforces import GBSAOBC2Force
from pdbfixer import PDBFixer

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
STRUCT = os.path.join(ROOT, "data", "structures")

ONE2THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}
THREE2ONE = {v: k for k, v in ONE2THREE.items()}
BACKBONE = ("N", "CA", "C", "O")
FF_FILES = ("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")
KCAL = 4.184


# --------------------------------------------------------------------------
# small geometry helpers
# --------------------------------------------------------------------------
def _dihedral(p0, p1, p2, p3):
    p0, p1, p2, p3 = (np.asarray(x, float) for x in (p0, p1, p2, p3))
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def _place(a, b, c, bond, angle_deg, dihedral_deg):
    """Atom d bonded to c with angle(b,c,d) and dihedral(a,b,c,d) (NeRF)."""
    a, b, c = (np.asarray(x, float) for x in (a, b, c))
    ang, dih = math.radians(angle_deg), math.radians(dihedral_deg)
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d2 = np.array([-bond * math.cos(ang), bond * math.sin(ang) * math.cos(dih),
                   bond * math.sin(ang) * math.sin(dih)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


def _atom_line(serial, name, resname, chain, resseq, xyz, element, record="ATOM  "):
    name_f = f" {name:<3}" if len(element) == 1 and len(name) < 4 else f"{name:<4}"
    return (f"{record}{serial:5d} {name_f} {resname:>3} {chain}{resseq:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2}\n")


# --------------------------------------------------------------------------
# template parsing
# --------------------------------------------------------------------------
def read_template(pdb_id="2fib", rec_chain="A", pep_chain="B"):
    """-> (receptor, ca_ion_xyz, peptide)

    receptor: dict resseq -> (resname, {atom: (xyz, element)}); altloc A/blank only
    peptide:  list of (resname, resseq, {atom: (xyz, element)})
    """
    rec, pep, ca = {}, {}, None
    for line in open(os.path.join(STRUCT, f"{pdb_id}.pdb")):
        if not line.startswith(("ATOM", "HETATM")) or line[16] not in (" ", "A"):
            continue
        resn, ch, name = line[17:20].strip(), line[21], line[12:16].strip()
        if resn == "HOH":
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        el = line[76:78].strip() or name[0]
        seq = int(line[22:26])
        if ch == rec_chain and resn == "CA" and line.startswith("HETATM"):
            ca = xyz
        elif ch == rec_chain:
            rec.setdefault(seq, (resn, {}))[1][name] = (xyz, el)
        elif ch == pep_chain:
            pep.setdefault(seq, (resn, {}))[1][name] = (xyz, el)
    peptide = [(pep[k][0], k, pep[k][1]) for k in sorted(pep)]
    return rec, ca, peptide


def _heavy_xyz(atoms):
    return np.array([x for n, (x, e) in atoms.items() if e != "H"])


def select_shell(rec, ca, peptide, shell):
    """Residue numbers within ``shell`` A of the template peptide (+Ca2+ site)."""
    pep_xyz = np.vstack([_heavy_xyz(a) for _, _, a in peptide])
    keep = set()
    for seq, (resn, atoms) in rec.items():
        d = np.linalg.norm(_heavy_xyz(atoms)[:, None, :] - pep_xyz[None], axis=2).min()
        if d <= shell:
            keep.add(seq)
    use_ca = False
    if ca is not None:
        ca_ligands = {seq for seq, (resn, atoms) in rec.items()
                      if np.linalg.norm(_heavy_xyz(atoms) - ca, axis=1).min() <= 3.0}
        if keep & ca_ligands:           # any ligand in the shell -> keep whole site
            keep |= ca_ligands
            use_ca = True
    ordered = sorted(keep)              # fill gaps of <= 2 residues
    for a, b in zip(ordered, ordered[1:]):
        if 1 < b - a <= 3:
            keep |= set(range(a + 1, b))
    return sorted(keep), use_ca


def segments(resids):
    segs, cur = [], [resids[0]]
    for r in resids[1:]:
        if r == cur[-1] + 1:
            cur.append(r)
        else:
            segs.append(cur)
            cur = [r]
    segs.append(cur)
    return segs


# --------------------------------------------------------------------------
# peptide construction
# --------------------------------------------------------------------------
def build_peptide_lines(sequence, peptide_t, start_serial, chain="B"):
    """Capped peptide as PDB lines.  Backbone (+CB where both residues have it)
    comes from the template residue at the same position; unchanged residues are
    copied whole.  Deletions drop trailing template residues; extensions are not
    supported.  PDBFixer adds the missing side-chain atoms afterwards."""
    assert len(sequence) <= len(peptide_t), "extensions are not supported"
    lines, serial = [], start_serial
    for i, aa in enumerate(sequence):
        resn_t, _, atoms_t = peptide_t[i]
        new = ONE2THREE[aa]
        keep = set(BACKBONE)
        if new == resn_t:
            keep |= set(atoms_t)
        elif new != "GLY" and "CB" in atoms_t:
            keep.add("CB")
        for name, (xyz, el) in atoms_t.items():
            if name == "OXT" or name not in keep:
                continue
            lines.append(_atom_line(serial, name, new, chain, i + 1, xyz, el))
            serial += 1
    last = peptide_t[len(sequence) - 1][2]
    n_prev, ca_, c_, o_ = (last[k][0] for k in ("N", "CA", "C", "O"))
    psi_o = _dihedral(n_prev, ca_, c_, o_)      # cap N is coplanar with C=O
    n_cap = _place(n_prev, ca_, c_, 1.335, 116.2, psi_o + 180.0)
    me_c = _place(ca_, c_, n_cap, 1.449, 121.9, 180.0)
    r = len(sequence) + 1
    lines.append(_atom_line(serial, "N", "NME", chain, r, n_cap, "N"))
    lines.append(_atom_line(serial + 1, "C", "NME", chain, r, me_c, "C"))
    return lines, serial + 2


class Model:
    def __init__(self, **kw):
        self.__dict__.update(kw)


def build_model(sequence="GPRP", shell=10.0, rec_mutations=None, ph=7.4, pdb_id="2fib"):
    """Build the pocket model.  ``rec_mutations``: {resseq: one-letter aa}
    applied to the receptor (side chain rebuilt by PDBFixer)."""
    rec, ca, peptide = read_template(pdb_id)
    resids, use_ca = select_shell(rec, ca, peptide, shell)
    rec_mutations = rec_mutations or {}
    for r in rec_mutations:
        assert r in resids, f"mutated residue {r} lies outside the {shell} A shell"
    lines, serial = [], 1
    for seg in segments(resids):
        a, b = seg[0], seg[-1]
        chain_lines = []
        if a - 1 in rec:                 # ACE from the real preceding residue
            for src, dst in (("CA", "CH3"), ("C", "C"), ("O", "O")):
                xyz, el = rec[a - 1][1][src]
                chain_lines.append(("ACE", a - 1, dst, xyz, el))
        for r in seg:
            resn, atoms = rec[r]
            new = ONE2THREE[rec_mutations[r]] if r in rec_mutations else resn
            keep = set(BACKBONE)
            if new == resn:
                keep |= set(atoms)
            elif new != "GLY" and "CB" in atoms:
                keep.add("CB")
            for name, (xyz, el) in atoms.items():
                if name in keep and name != "OXT":
                    chain_lines.append((new, r, name, xyz, el))
        if b + 1 in rec:                 # NME from the real following residue
            for src, dst in (("N", "N"), ("CA", "C")):
                xyz, el = rec[b + 1][1][src]
                chain_lines.append(("NME", b + 1, dst, xyz, el))
        for resn, r, name, xyz, el in chain_lines:
            lines.append(_atom_line(serial, name, resn, "A", r, xyz, el))
            serial += 1
        lines.append("TER\n")
    if use_ca:
        lines.append(_atom_line(serial, "CA", "CA", "C", 1, ca, "Ca", record="HETATM"))
        serial += 1
        lines.append("TER\n")
    pep_lines, serial = build_peptide_lines(sequence, peptide, serial)
    lines += pep_lines + ["TER\n", "END\n"]
    tmp = os.path.join(ROOT, "data", f"_tmp_build_{os.getpid()}.pdb")
    with open(tmp, "w") as fh:
        fh.writelines(lines)
    fx = PDBFixer(filename=tmp)
    fx.findMissingResidues()
    fx.missingResidues = {}
    fx.findNonstandardResidues()
    fx.replaceNonstandardResidues()
    fx.findMissingAtoms()
    fx.addMissingAtoms()
    fx.addMissingHydrogens(ph)
    os.remove(tmp)
    top, pos = fx.topology, fx.positions
    pep_chain = list(top.chains())[-1]
    pep_idx = [a.index for r in pep_chain.residues() for a in r.atoms()]
    n_rec = min(pep_idx)
    assert pep_idx == list(range(n_rec, top.getNumAtoms())), "peptide must be last"
    atoms = list(top.atoms())
    anchors = sorted({a.index for a in atoms[:n_rec]
                      if (a.name == "CA" and a.residue.name not in ("ACE", "NME"))
                      or (a.residue.name == "ACE" and a.name == "CH3")
                      or (a.residue.name == "NME" and a.name == "C")})
    pep_bb = [a.index for a in atoms[n_rec:]
              if a.name in ("N", "CA", "C") and a.residue.name != "NME"]
    return Model(topology=top, positions=pos, n_rec=n_rec, anchors=anchors,
                 pep_idx=pep_idx, pep_bb=pep_bb, sequence=sequence, shell=shell,
                 resids=resids, use_ca=use_ca, rec_mutations=rec_mutations)


# --------------------------------------------------------------------------
# force field / simulation
# --------------------------------------------------------------------------
def make_system(topology, hmr=True):
    ff = app.ForceField(*FF_FILES)
    kw = dict(nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    if hmr:
        kw["hydrogenMass"] = 3.0 * unit.amu
    system = ff.createSystem(topology, **kw)
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    nb.setReactionFieldDielectric(1.0)
    gb = mm.GBSAOBCForce()
    gb.setNonbondedMethod(mm.GBSAOBCForce.NoCutoff)
    gb.setSoluteDielectric(1.0)
    gb.setSolventDielectric(78.5)
    gb.setForceGroup(2)             # group 0 = gas phase, 1 = restraints, 2 = GB+SA
    for i, (r, s) in enumerate(GBSAOBC2Force.getStandardParameters(topology)):
        gb.addParticle(nb.getParticleParameters(i)[0], r, s)
    system.addForce(gb)
    return system


def add_anchor_restraints(system, positions, atom_ids, k_kcal_a2):
    f = mm.CustomExternalForce("0.5*kr*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    f.addGlobalParameter("kr", k_kcal_a2 * KCAL * 100.0)
    for p in ("x0", "y0", "z0"):
        f.addPerParticleParameter(p)
    for i in atom_ids:
        f.addParticle(i, positions[i].value_in_unit(unit.nanometer))
    f.setForceGroup(1)
    system.addForce(f)


def add_leash(system, positions, pep_bb, radius_a=6.0, k_kcal_a2=10.0):
    """Flat-bottom restraint on the peptide-backbone centroid (no force inside)."""
    xyz = np.array(positions.value_in_unit(unit.nanometer))
    c0 = xyz[pep_bb].mean(axis=0)
    f = mm.CustomCentroidBondForce(
        1, "0.5*kl*step(dd-r0)*(dd-r0)^2; dd=sqrt((x1-cx)^2+(y1-cy)^2+(z1-cz)^2)")
    f.addGlobalParameter("kl", k_kcal_a2 * KCAL * 100.0)
    f.addGlobalParameter("r0", radius_a / 10.0)
    for name, v in zip(("cx", "cy", "cz"), c0):
        f.addGlobalParameter(name, float(v))
    f.addBond([f.addGroup(list(pep_bb))], [])
    f.setForceGroup(1)
    system.addForce(f)


def run_md(model, seed, equil_ps=80, prod_ps=200, save_ps=5, threads=1,
           temperature=300.0, k_rest=1.0, anneal=True, anneal_ps=5, leash_a=6.0):
    """Minimise -> anneal (450 K, randomises peptide side chains) -> equilibrate
    -> production.  Returns (frames in nm, meta)."""
    system = make_system(model.topology)
    add_anchor_restraints(system, model.positions, model.anchors, k_rest)
    add_leash(system, model.positions, model.pep_bb, leash_a)
    integ = mm.LangevinMiddleIntegrator(temperature * unit.kelvin, 1.0 / unit.picosecond,
                                        4.0 * unit.femtosecond)
    integ.setRandomNumberSeed(seed)
    ctx = mm.Context(system, integ, mm.Platform.getPlatformByName("CPU"),
                     {"Threads": str(threads)})
    ctx.setPositions(model.positions)
    mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 1000)
    spp = 250                                     # steps per ps at 4 fs
    t0 = time.time()
    ctx.setVelocitiesToTemperature(100 * unit.kelvin, seed)
    ramp = (150, 250, 350, 450, 400, 350, 300) if anneal else (150, 200, 250, 300)
    for T in ramp:
        integ.setTemperature(T * unit.kelvin)
        integ.step(int(anneal_ps * spp))
    integ.setTemperature(temperature * unit.kelvin)
    integ.step(int(equil_ps * spp))
    n = int(prod_ps / save_ps)
    frames = np.empty((n, model.topology.getNumAtoms(), 3), dtype=np.float32)
    for i in range(n):
        integ.step(int(save_ps * spp))
        frames[i] = ctx.getState(getPositions=True).getPositions(
            asNumpy=True).value_in_unit(unit.nanometer)
    meta = dict(seed=seed, equil_ps=equil_ps, prod_ps=prod_ps, anneal=anneal,
                wall_s=round(time.time() - t0, 1))
    return frames, meta


# --------------------------------------------------------------------------
# single-trajectory MM/GBSA (same OBC-II model as the dynamics)
# --------------------------------------------------------------------------
class GBSAScorer:
    def __init__(self, model):
        top, n = model.topology, model.topology.getNumAtoms()
        self.sets = {"cx": list(range(n)), "rc": list(range(model.n_rec)),
                     "lg": list(range(model.n_rec, n))}
        self.ctx = {}
        for key, ids in self.sets.items():
            mod = app.Modeller(top, model.positions)
            keep = set(ids)
            mod.delete([a for a in mod.topology.atoms() if a.index not in keep])
            self.ctx[key] = mm.Context(make_system(mod.topology, hmr=False),
                                       mm.VerletIntegrator(0.001),
                                       mm.Platform.getPlatformByName("CPU"), {"Threads": "1"})

    def energies(self, frame_nm):
        """-> {set: (E_gas, E_gb)} in kcal/mol (restraint group excluded)."""
        out = {}
        for key, ids in self.sets.items():
            c = self.ctx[key]
            c.setPositions(frame_nm[ids])
            g = c.getState(getEnergy=True, groups={0}).getPotentialEnergy()
            w = c.getState(getEnergy=True, groups={2}).getPotentialEnergy()
            out[key] = (g.value_in_unit(unit.kilocalories_per_mole),
                        w.value_in_unit(unit.kilocalories_per_mole))
        return out

    def score(self, frames):
        """-> dict of per-frame arrays: dg (total), dg_gas, dg_gb."""
        e = [self.energies(f) for f in frames]
        arr = {k: np.array([x[k] for x in e]) for k in ("cx", "rc", "lg")}   # (n,2)
        gas = arr["cx"][:, 0] - arr["rc"][:, 0] - arr["lg"][:, 0]
        gb = arr["cx"][:, 1] - arr["rc"][:, 1] - arr["lg"][:, 1]
        return dict(dg=gas + gb, dg_gas=gas, dg_gb=gb)


def write_pdb(model, frame_nm, path):
    with open(path, "w") as fh:
        app.PDBFile.writeFile(model.topology, np.asarray(frame_nm) * unit.nanometer,
                              fh, keepIds=True)
