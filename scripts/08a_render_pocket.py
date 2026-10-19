"""Render the hole-a pocket of 2FIB with the bound GPRP (Figure 1A) using PyMOL.

Ubuntu's PyMOL 2.5 targets Python 3.12 and still imports the removed `imp`
module, so run with python3.12 and a tiny `imp` shim on PYTHONPATH, e.g.
  PYTHONPATH=<dir containing imp.py> /usr/bin/python3.12 scripts/08a_render_pocket.py
(the shim only needs `find_module`; see docs/ for the 8-line file).
"""
import os
import sys

sys.path.insert(0, "/usr/lib/python3/dist-packages")
import pymol2  # noqa: E402

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_p = pymol2.PyMOL()
_p.start()
cmd = _p.cmd
cmd.load(os.path.join(root, "data/structures/2fib.pdb"), "m")
cmd.remove("solvent")
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("ray_shadows", 0)
cmd.set("depth_cue", 0)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("stick_radius", 0.17)
cmd.set("label_size", 22)
cmd.set("label_font_id", 7)
cmd.set("label_color", "black")
cmd.set("dash_gap", 0.25)
cmd.set("dash_radius", 0.07)
cmd.set("dash_color", "black")

cmd.select("pep", "chain B")
cmd.select("rec", "chain A and polymer")
cmd.select("ion", "resn CA and chain A and hetatm")
cmd.select("pocket", "(byres rec within 4.2 of pep)")
cmd.hide("everything")
cmd.show("cartoon", "rec within 13 of pep")
cmd.set("cartoon_transparency", 0.55)
cmd.set("cartoon_loop_radius", 0.12)
cmd.color("grey80", "rec")
cmd.show("sticks", "pocket")
cmd.set("stick_radius", 0.11, "pocket")
cmd.color("lightblue", "pocket and elem C")
cmd.show("sticks", "pep")
cmd.set("stick_radius", 0.27, "pep")
cmd.color("orange", "pep and elem C")
for _e, _c in (("N", "blue"), ("O", "red"), ("S", "yellow")):
    cmd.color(_c, f"elem {_e}")
cmd.show("spheres", "ion")
cmd.set("sphere_scale", 0.5)
cmd.color("forest", "ion")
cmd.distance("hb", "pep", "pocket", 3.4, mode=2)
cmd.hide("labels", "hb")
cmd.set("dash_radius", 0.06)

# colour-code key residues (legend is provided by the contact map in Figure 1B)
cmd.color("salmon", "pocket and resi 364+330+297+301 and elem C")     # acidic
cmd.color("lightteal", "pocket and resi 329+338+340+375 and elem C")  # polar / basic
cmd.color("sand", "pocket and resi 295+322+305 and elem C")           # aromatic / Thr
cmd.set("stick_radius", 0.16, "pocket and resi 364+330+329+297+301+338+340+375+295+322+305")
cmd.orient("pep or pocket")
cmd.zoom("pep or pocket", 2)
cmd.turn("y", 35)
cmd.turn("x", -25)
cmd.png(os.path.join(root, "figures/src/fig1a_pocket.png"), width=1600, height=1300, dpi=300, ray=1)
_p.stop()
