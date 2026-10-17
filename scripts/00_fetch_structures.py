"""Fetch candidate fibrinogen/fibrin fragment-D structures from the wwPDB
AWS Open Data mirror and report title / ligand / chain content.

Only a small set of mirror hosts is reachable from the analysis sandbox
(rcsb.org itself is not), so every structure used in the paper is pulled from
the public wwPDB snapshot bucket below.  The snapshot date is recorded.
"""
import gzip, io, sys, urllib.request, pathlib, json

SNAPSHOT = "20240101"
BASE = ("https://pdbsnapshots.s3.us-west-2.amazonaws.com/%s/pub/pdb/data/"
        "structures/divided/pdb/%s/pdb%s.ent.gz")
OUT = pathlib.Path(__file__).resolve().parent.parent / "data" / "structures"


def fetch(pdb_id: str) -> pathlib.Path:
    pdb_id = pdb_id.lower()
    dest = OUT / f"{pdb_id}.pdb"
    if dest.exists():
        return dest
    url = BASE % (SNAPSHOT, pdb_id[1:3], pdb_id)
    with urllib.request.urlopen(url, timeout=60) as r:
        dest.write_bytes(gzip.decompress(r.read()))
    return dest


def summarize(path: pathlib.Path) -> dict:
    title, compnd, het, remark2, seqres = [], [], set(), "", {}
    for line in path.read_text().splitlines():
        rec = line[:6].strip()
        if rec == "TITLE":
            title.append(line[10:].strip())
        elif rec == "COMPND":
            compnd.append(line[10:].strip())
        elif rec == "HET":
            het.add(line[7:10].strip())
        elif rec == "REMARK" and line[7:10] == "2 " and "RESOLUTION" in line:
            remark2 = line[11:].strip()
        elif rec == "SEQRES":
            seqres.setdefault(line[11], []).extend(line[19:].split())
    short = {c: "".join(s) if len(s) < 40 else f"<{len(s)} res>"
             for c, s in seqres.items()}
    return dict(id=path.stem.upper(), title=" ".join(title), res=remark2,
                het=sorted(het), chains={c: ("-".join(s) if len(s) <= 12 else f"{len(s)} res")
                                         for c, s in seqres.items()})


if __name__ == "__main__":
    ids = sys.argv[1:]
    for i in ids:
        try:
            s = summarize(fetch(i))
            print(json.dumps(s))
        except Exception as e:
            print(json.dumps(dict(id=i.upper(), error=str(e))))
