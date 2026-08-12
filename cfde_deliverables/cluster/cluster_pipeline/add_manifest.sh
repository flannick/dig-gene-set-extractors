#!/usr/bin/env bash
# Add/refresh MANIFEST.md + manifest.json for a completed gene-set submission. No re-run.
# GENERIC contract-format introspection. Accepts EITHER:
#   add_manifest.sh <submission.zip>        # unzip -> add manifest -> re-zip in place -> refresh .md5
#   add_manifest.sh <geneset_dir> <name>    # write manifest into dir -> zip to SDG as <name>.zip
# Auto-skips bundles with no gene sets (e.g. script zips). Pure python3 stdlib + zip/unzip.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

gen_manifest() {  # $1 = dir to scan & write into, $2 = display name ; exit 3 if no gene sets
  python3 - "$1" "$2" <<'PY'
import sys, os, glob, json, datetime, statistics
DIR, NAME = sys.argv[1], sys.argv[2]
metas=[]
for f in glob.glob(os.path.join(DIR,"**","geneset.meta.json"), recursive=True):
    try: metas.append((f,json.load(open(f))))
    except Exception: pass
if not metas:
    print("NO_GENESETS"); sys.exit(3)
libs={}; orgs=set(); sources=set(); fundings=set(); pub=set(); ng=[]; methods=set(); ks=set()
for f,m in metas:
    libs[m.get("library","?")]=libs.get(m.get("library","?"),0)+1
    if m.get("organism"): orgs.add(str(m["organism"]))
    if m.get("source"): sources.add(m["source"])
    if isinstance(m.get("n_genes"),int): ng.append(m["n_genes"])
    if m.get("method"): methods.add(str(m["method"]))
    if m.get("K") is not None: ks.add(str(m["K"]))
    pv=os.path.join(os.path.dirname(f),"geneset.provenance.json")
    if os.path.exists(pv):
        try:
            p=json.load(open(pv))
            if p.get("funding"): fundings.add(p["funding"])
            pub.add(bool(p.get("public")))
        except Exception: pass
n=len(metas); gmin,gmed,gmax=(min(ng),int(statistics.median(ng)),max(ng)) if ng else (0,0,0)
obj={"submission":NAME,"generated":datetime.date.today().isoformat(),"n_genesets":n,
     "genes_per_set":{"min":gmin,"median":gmed,"max":gmax},"libraries":libs,
     "organisms":sorted(orgs),"methods":sorted(methods),"K":sorted(ks),
     "publicly_accessible":(pub=={True}),"funding":sorted(fundings),"sources":sorted(sources)}
json.dump(obj, open(os.path.join(DIR,"manifest.json"),"w"), indent=1)
libtbl="\n".join(f"- {k}: {v}" for k,v in sorted(libs.items(),key=lambda x:-x[1]))
md=f"""# Submission MANIFEST — {NAME}

Generated: {obj['generated']} (auto-generated from the completed gene sets)
Gene sets: {n} | genes per set: min {gmin} / median {gmed} / max {gmax}

## Libraries
{libtbl}

## Organisms
{", ".join(sorted(orgs)) or "n/a"}
"""
if methods:
    md+=f"\n## Method\n- {', '.join(sorted(methods))}"+(f" (K={', '.join(sorted(ks))})" if ks else "")+"\n"
md+="\n## Compliance\n- Publicly accessible: "+("YES (all sets flagged public)" if pub=={True} else "MIXED/UNKNOWN — check per-set provenance")+"\n- Funding:\n"
md+=("\n".join(f"- {x}" for x in sorted(fundings)) or "- (see per-set provenance)")
md+="\n\n## Sources (deduplicated; up to 25)\n"+("\n".join(f"- {s}" for s in sorted(sources)[:25]) or "- (see per-set provenance)")
md+="\n\n## Format\nEach set: geneset.tsv / genesets.gmt / geneset.meta.json / geneset.provenance.json.\nMachine-readable summary: manifest.json (this folder).\n"
open(os.path.join(DIR,"MANIFEST.md"),"w").write(md)
print(f"OK {n} sets, {len(libs)} libraries")
PY
}

T=${1:?usage: add_manifest.sh <submission.zip>  OR  <geneset_dir> <name>}
if [[ "$T" == *.zip ]]; then
  [ -f "$T" ] || { echo "no such zip: $T"; exit 1; }
  NAME=$(basename "$T" .zip); TMP=$(mktemp -d)
  unzip -q "$T" -d "$TMP"
  if gen_manifest "$TMP" "$NAME"; then
    NEW="$TMP.repack.zip"; ( cd "$TMP" && zip -rq "$NEW" . ) && mv "$NEW" "$T"
    md5sum "$T" | tee "$T.md5"
    echo "UPDATED $T (added MANIFEST.md + manifest.json)"
  else
    echo "SKIP $T — no gene sets inside (probably a scripts bundle); left unchanged"
  fi
  rm -rf "$TMP"
else
  DIR="$T"; NAME=${2:?need <name> when first arg is a directory}; [ -d "$DIR" ] || { echo "no such dir: $DIR"; exit 1; }
  if gen_manifest "$DIR" "$NAME"; then
    bash "$HERE/package_submission.sh" "$DIR" "$NAME"; echo "re-zipped $NAME with MANIFEST.md + manifest.json"
  else
    echo "SKIP $DIR — no gene sets found"
  fi
fi
