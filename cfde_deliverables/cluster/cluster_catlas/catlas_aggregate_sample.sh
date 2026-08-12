#!/bin/bash
# Aggregate ONE CATLAS fragment file -> per-cell-type promoter fragment counts.
# Args: $1 = fragment .bed.gz path.  Env: WORK (has promoter_bins.tsv + metadata), OUTDIR (partials).
# Sample-name normalization: GSM#####_<sample>_repN_fragments.bed.gz  ->  metadata sample "<sample>_N".
set -euo pipefail
FRAG="$1"
WORK="${WORK:?set WORK}"; OUTDIR="${OUTDIR:?set OUTDIR}"
META="${META:-$WORK/GSE184462_metadata.tsv.gz}"
BINSZ="${BIN_BP:-1000}"
mkdir -p "$OUTDIR"
base=$(basename "$FRAG" _fragments.bed.gz)          # GSM#####_adipose_omentum_SM-ADYHB_rep1
nogsm=${base#GSM*_}                                  # adipose_omentum_SM-ADYHB_rep1
sample=$(echo "$nogsm" | sed -E 's/_rep([0-9]+)$/_\1/')   # adipose_omentum_SM-ADYHB_1
out="$OUTDIR/${sample}.partial.tsv"
tmpmap="$OUTDIR/${sample}.map.tsv"

# barcode -> cell type map for THIS sample (col1 cellID = "sample+barcode"; col7 = cell type)
gzip -dc "$META" | awk -F'\t' -v s="$sample" 'NR>1 && $2==s { n=split($1,a,"+"); print a[n]"\t"$7 }' > "$tmpmap"
ncells=$(wc -l < "$tmpmap")
echo "[$sample] barcodes in metadata: $ncells" >&2
if [ "$ncells" -eq 0 ]; then echo "[$sample] WARN no metadata match, skipping" >&2; : > "$out"; rm -f "$tmpmap"; exit 0; fi

# stream: map(FS \t) + bins(FS \t) + fragments(FS \t) ; count 1 per fragment overlapping a promoter bin
awk -F'\t' -v BINSZ="$BINSZ" '
  FNR==1{fno++}
  fno==1{ CT[$1]=$2; next }                                   # barcode -> cell type
  fno==2{ BING[$1]=$2; next }                                 # chrom:bin -> gene1;gene2
  { ct=CT[$4]; if(ct=="") next;
    b1=int($2/BINSZ); b2=int($3/BINSZ);
    for(b=b1;b<=b2;b++){ key=$1":"b;
      if(key in BING){ m=split(BING[key],gg,";"); for(i=1;i<=m;i++) CNT[ct SUBSEP gg[i]]++ } }
  }
  END{ for(k in CNT){ split(k,p,SUBSEP); print p[1]"\t"p[2]"\t"CNT[k] } }
' "$tmpmap" "$WORK/promoter_bins.tsv" <(gzip -dc "$FRAG") > "$out"
rm -f "$tmpmap"
echo "[$sample] wrote $(wc -l < "$out") celltype-gene rows -> $out" >&2
