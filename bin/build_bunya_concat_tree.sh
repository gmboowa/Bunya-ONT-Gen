#!/usr/bin/env bash
# build_bunya_concat_tree.sh
# Single phylogeny from three Bunyamwera segments (concatenated).
# - Per-sample consensus per segment (headers: SAMPLE|SEGMENT)
# - MAFFT per-segment alignments
# - Robust concatenation by SAMPLE NAME (not record order)
# - IQ-TREE (UFboot + SH-aLRT) or FastTree (-boot) fallback
#
# Install (conda/mamba):
#   mamba install -c bioconda -c conda-forge bcftools samtools mafft seqkit iqtree fasttree entrez-direct
#
# Usage:
#   ./build_bunya_concat_tree.sh -l vcfs.txt -o outdir [--auto-ref] [--segments A,B,C] [--threads N] [--boot 1000] [--iupac]
#   vcfs.txt: one VCF per line OR "sample<TAB>/path/to/sample.vcf[.gz]"
set -euo pipefail
trap 'echo "ERROR: Command failed: ${BASH_COMMAND} (line $LINENO)" >&2' ERR

# ---------- defaults ----------
LIST=""; OUTDIR=""
THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)"
BOOT=1000
AUTO_REF=0
SEGMENTS=""
REFL=""; REFM=""; REFS=""
IUPAC=0
MASK_BED=""

die(){ echo "ERROR: $*" >&2; exit 1; }

usage(){ cat <<EOF
Usage:
  $0 -l <vcf_list.txt> -o <outdir> [options]

Required:
  -l, --list FILE          VCF list: one per line OR "sample<TAB>path"
  -o, --outdir DIR         Output directory

References (choose one):
  --auto-ref               Fetch FASTA via NCBI efetch for each detected segment
  --refL FILE --refM FILE --refS FILE
                           Provide L/M/S FASTAs; headers rewritten to match VCF contig IDs

Segments:
  --segments A,B,C         Three segment IDs as they appear in the VCFs (overrides auto-detect)

Other:
  --threads N              Threads for MAFFT/IQ-TREE (default: auto)
  --boot N                 Bootstrap replicates (IQ-TREE UFboot & SH-aLRT; FastTree -boot) [default: 1000]
  --mask-bed BED           BED of positions to mask (N) before consensus
  --iupac                  Keep IUPAC ambiguity codes (else pick haplotype -H 1)
EOF
}

# ---------- parse args ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -l|--list) LIST="$2"; shift 2;;
    -o|--outdir) OUTDIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --boot) BOOT="$2"; shift 2;;
    --auto-ref) AUTO_REF=1; shift;;
    --segments) SEGMENTS="$2"; shift 2;;
    --refL) REFL="$2"; shift 2;;
    --refM) REFM="$2"; shift 2;;
    --refS) REFS="$2"; shift 2;;
    --mask-bed) MASK_BED="$2"; shift 2;;
    --iupac) IUPAC=1; shift;;
    -h|--help) usage; exit 0;;
    *) die "Unknown argument: $1 (see -h)";;
  esac
done

[[ -z "$LIST" || -z "$OUTDIR" ]] && { usage; exit 1; }
[[ -f "$LIST" ]] || die "List file not found: $LIST"
mkdir -p "$OUTDIR"/{refs,tmp,cons,aln,trees,logs}

# ---------- deps ----------
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }
for tool in bcftools bgzip tabix samtools mafft seqkit awk grep sed cut sort uniq; do need "$tool"; done
[[ "$AUTO_REF" -eq 0 ]] || need efetch

IQPROG=""; IQVER=0
if command -v iqtree2 >/dev/null 2>&1; then IQPROG="$(command -v iqtree2)"; IQVER=2;
elif command -v iqtree >/dev/null 2>&1; then IQPROG="$(command -v iqtree)"; IQVER=1;
fi

FASTTREE=""
if command -v fasttreeMP >/dev/null 2>&1; then FASTTREE="$(command -v fasttreeMP)";
elif command -v FastTreeMP >/dev/null 2>&1; then FASTTREE="$(command -v FastTreeMP)";
elif command -v fasttree >/dev/null 2>&1; then FASTTREE="$(command -v fasttree)";
elif command -v FastTree >/dev/null 2>&1; then FASTTREE="$(command -v FastTree)";
fi

if [[ -z "$IQPROG" && -z "$FASTTREE" ]]; then
  die "Neither IQ-TREE (iqtree2/iqtree) nor FastTree found. Install: mamba install -c bioconda iqtree fasttree"
fi

THREAD_FLAG=""; SAFE_THREAD_FLAG=""
MODEL_FLAG=""; UFBOOT_FLAG=""; PART_FLAG="-p"
if [[ -n "$IQPROG" ]]; then
  if "$IQPROG" -h 2>&1 | grep -q '\-T'; then THREAD_FLAG="-T $THREADS"; SAFE_THREAD_FLAG="-T 1"; else THREAD_FLAG="-nt $THREADS"; SAFE_THREAD_FLAG="-nt 1"; fi
  if [[ $IQVER -eq 2 ]]; then UFBOOT_FLAG="-B $BOOT -alrt $BOOT"; MODEL_FLAG="MFP+MERGE"; PART_FLAG="-p";
  else UFBOOT_FLAG="-bb $BOOT -alrt $BOOT"; MODEL_FLAG="TEST"; PART_FLAG=$([[ "$("$IQPROG" -h 2>&1 | grep -c '\-p ')" -gt 0 ]] && echo "-p" || echo "-spp"); fi
fi

# sed -i helper (GNU vs BSD)
sed_inplace(){ if sed --version >/dev/null 2>&1; then sed -E -i "$1" "$2"; else sed -E -i '' "$1" "$2"; fi; }

# ---------- read sample list ----------
declare -a SAMPLES VCFPATHS
while IFS= read -r line; do
  [[ -z "$line" || "$line" =~ ^# ]] && continue
  if echo "$line" | grep -q $'\t'; then
    sid="$(echo "$line" | cut -f1)"; vcf="$(echo "$line" | cut -f2-)"
  else
    vcf="$line"; sid="$(basename "$vcf")"; sid="${sid%.*}"
  fi
  [[ -f "$vcf" ]] || die "VCF not found: $vcf"
  SAMPLES+=("$sid"); VCFPATHS+=("$vcf")
done < "$LIST"
[[ ${#SAMPLES[@]} -ge 2 ]] || die "Need at least 2 VCFs."

# ---------- detect segments ----------
IFS=',' read -r -a SEG_ARR <<< "${SEGMENTS:-}"
if [[ ${#SEG_ARR[@]} -ne 3 ]]; then
  echo "Auto-detecting segment IDs from first VCF..."
  VCF0="${VCFPATHS[0]}"; LC_ALL=C
  header_ids_str="$(grep -E '^##contig=<ID=' "$VCF0" | sed -E 's/.*ID=([^,>]+).*/\1/')" || true
  header_count="$(printf "%s\n" "${header_ids_str:-}" | grep -c . || true)"
  if [[ "${header_count:-0}" -lt 3 ]]; then
    echo "Falling back to CHROM frequencies..."
    header_ids_str="$(grep -v '^#' "$VCF0" | cut -f1 | sort | uniq -c | sort -nr | awk '{print $2}' | head -n 3)" || true
    header_count="$(printf "%s\n" "${header_ids_str:-}" | grep -c . || true)"
  fi
  [[ "${header_count:-0}" -lt 3 ]] && die "Could not detect 3 segment IDs"
  SEG_ARR=(); while IFS= read -r id; do [[ -n "$id" ]] && SEG_ARR+=("$id"); done <<EOF
${header_ids_str:-}
EOF
  SEG_ARR=("${SEG_ARR[@]:0:3}")
fi
echo "Using segments: ${SEG_ARR[*]}"

# ---------- references ----------
for seg in "${SEG_ARR[@]}"; do
  OUTFA="$OUTDIR/refs/${seg}.fa"
  if [[ "$AUTO_REF" -eq 1 ]]; then
    echo "Fetching reference for $seg via efetch..."
    efetch -db nucleotide -format fasta -id "$seg" > "$OUTFA"
    grep -q '^>' "$OUTFA" || die "No FASTA fetched for $seg"
    sed_inplace "1s/^>.*/>${seg}/" "$OUTFA"
  else
    [[ -n "$REFL" && -n "$REFM" && -n "$REFS" ]] || die "Provide --auto-ref OR all of --refL/--refM/--refS"
    case "$seg" in
      "${SEG_ARR[0]}") SRC="$REFL";;
      "${SEG_ARR[1]}") SRC="$REFM";;
      "${SEG_ARR[2]}") SRC="$REFS";;
      *) die "Internal segment mapping error";;
    esac
    [[ -f "$SRC" ]] || die "Reference FASTA not found: $SRC"
    awk -v id="$seg" 'BEGIN{p=0} /^>/{if(!p){print ">" id; p=1}else{print}} !/^>/{print}' "$SRC" > "$OUTFA"
  fi
  samtools faidx "$OUTFA"
done

# ---------- compress/index VCFs & sample names ----------
declare -a SRC_VCF_GZ SRC_SAMPLES
for i in "${!SAMPLES[@]}"; do
  sid="${SAMPLES[$i]}"; vcf="${VCFPATHS[$i]}"
  if [[ "$vcf" =~ \.vcf\.gz$ ]]; then
    gz="$vcf"; [[ -f "${gz}.tbi" ]] || tabix -p vcf "$gz"
  else
    gz="$OUTDIR/tmp/${sid}.full.vcf.gz"
    bgzip -c "$vcf" > "$gz"; tabix -p vcf "$gz"
  fi
  SRC_VCF_GZ[$i]="$gz"
  sname="$(bcftools query -l "$gz" | head -n1 || true)"
  SRC_SAMPLES[$i]="$sname"
done

# ---------- helpers ----------
fasta_has_seq(){
  local fa="$1"
  [[ -s "$fa" ]] || return 1
  local n; n=$(awk 'BEGIN{c=0} /^>/ {next} {gsub(/[[:space:]]/,""); c+=length($0)} END{print c}' "$fa")
  [[ "${n:-0}" -ge 1 ]]
}
alignment_varsites(){
  local aln="$1"
  awk '
    BEGIN{FS=""}
    /^>/ {next}
    { for(i=1;i<=NF;i++){ c=toupper($i); if(c!~/[N\-]/){ seen[i SUBSEP c]=1 } } }
    END{
      for(k in seen){ split(k,a,SUBSEP); col=a[1]; cnt[col]++ }
      vs=0; for(col in cnt){ if(cnt[col]>1) vs++ } print vs
    }' "$aln"
}
alignment_nseq(){ grep -c '^>' "$1" 2>/dev/null || echo 0; }
alignment_unique(){ seqkit seq -w 0 "$1" 2>/dev/null | awk '/^>/{next}{print}' | sort | uniq | wc -l | tr -d ' ' || echo 0; }
sanitize_name(){ echo "$1" | sed -E 's/[[:space:]:,;()]+/_/g'; }
write_star_tree(){
  local aln="$1" out="$2"
  local names; names=$(grep '^>' "$aln" || true)
  [[ -n "$names" ]] || { echo "(X:0.0);" > "$out"; return 0; }
  local first=1 buf="("
  while IFS= read -r h; do
    [[ -z "$h" ]] && continue
    local n="${h#>}"; n="$(sanitize_name "$n")"
    if [[ $first -eq 1 ]]; then buf+="${n}:0.0"; first=0; else buf+=",${n}:0.0"; fi
  done <<< "$names"
  buf+=")"; echo "$buf;" > "$out"
}
rename_first_header(){ if sed --version >/dev/null 2>&1; then sed -i "1s/^>.*/>$2/" "$1"; else sed -i '' "1s/^>.*/>$2/" "$1"; fi; }
ensure_unique_headers(){
  local fa="$1"
  awk '
    /^>/ {h=$0; sub(/^>/,"",h); c[h]++; if(c[h]>1){printf(">%s_%d\n", h, c[h]);} else {print ">"h} next}
    {print}
  ' "$fa" > "${fa}.uniq.tmp"
  mv "${fa}.uniq.tmp" "$fa"
}

# ---------- normalize -> consensus ----------
make_consensus(){
  local seg="$1" sid="$2" gz="$3" sname="$4" reffa="$5" outfa="$6"
  local logv="$OUTDIR/logs/${sid}.${seg}.view.log"
  local logn="$OUTDIR/logs/${sid}.${seg}.norm.log"
  local logc="$OUTDIR/logs/${sid}.${seg}.cons.log"
  local tmpbcf="$OUTDIR/tmp/${sid}.${seg}.norm.bcf"

  if ! bcftools view -r "$seg" -Ou "$gz" 2> "$logv" \
    | bcftools view -v snps,indels -e 'ALT ~ "<"' -Ou \
    | bcftools norm -f "$reffa" --check-ref w -m -both -Ob -o "$tmpbcf" 2> "$logn"; then
    echo "WARN: extract/normalize failed for $sid/$seg; using reference" | tee -a "$logv" "$logn"
    cp "$reffa" "$outfa"
  else
    bcftools index -f "$tmpbcf" >/dev/null 2>>"$logn" || true
    CONS_ARGS=()
    [[ -n "$MASK_BED" ]] && CONS_ARGS+=( -m "$MASK_BED" )
    if [[ -n "$sname" ]]; then
      CONS_ARGS+=( -s "$sname" )
      if [[ "$IUPAC" -eq 1 ]]; then CONS_ARGS+=( -I ); else CONS_ARGS+=( -H 1 ); fi
    fi
    if ! bcftools consensus -f "$reffa" "${CONS_ARGS[@]}" "$tmpbcf" > "$outfa" 2> "$logc"; then
      echo "WARN: consensus failed for $sid/$seg; using reference" | tee -a "$logc"
      cp "$reffa" "$outfa"
    fi
  fi
  if ! fasta_has_seq "$outfa"; then
    echo "WARN: empty consensus for $sid/$seg; replacing with reference" | tee -a "$logc"
    cp "$reffa" "$outfa"
  fi
  rename_first_header "$outfa" "${sid}|${seg}"
}

# ---------- build consensus ----------
echo "==> Building consensus sequences ..."
CONS_DIR="$OUTDIR/cons"
for i in "${!SAMPLES[@]}"; do
  sid="${SAMPLES[$i]}"; gz="${SRC_VCF_GZ[$i]}"; sname="${SRC_SAMPLES[$i]}"
  echo "Processing $sid ..."
  for seg in "${SEG_ARR[@]}"; do
    make_consensus "$seg" "$sid" "$gz" "$sname" "$OUTDIR/refs/${seg}.fa" "$CONS_DIR/${sid}.${seg}.fa"
  done
done

# Ensure files exist
for sid in "${SAMPLES[@]}"; do
  for seg in "${SEG_ARR[@]}"; do
    fa="$CONS_DIR/${sid}.${seg}.fa"
    if ! fasta_has_seq "$fa"; then
      echo "WARN: repairing missing/empty $fa from reference"
      cp "$OUTDIR/refs/${seg}.fa" "$fa"; rename_first_header "$fa" "${sid}|${seg}"
    fi
  done
done

# ---------- per-segment alignments ----------
echo "==> Aligning per-segment consensuses ..."
ALN_DIR="$OUTDIR/aln"
for seg in "${SEG_ARR[@]}"; do
  : > "$ALN_DIR/${seg}.all.fa"
  for sid in "${SAMPLES[@]}"; do cat "$CONS_DIR/${sid}.${seg}.fa" >> "$ALN_DIR/${seg}.all.fa"; done
  ensure_unique_headers "$ALN_DIR/${seg}.all.fa"
  nseq=$(alignment_nseq "$ALN_DIR/${seg}.all.fa")
  [[ "$nseq" -ge 2 ]] || die "Not enough sequences to align for segment $seg (found $nseq)."
  mafft --thread "$THREADS" --auto "$ALN_DIR/${seg}.all.fa" > "$ALN_DIR/${seg}.aln.fa" 2> "$OUTDIR/logs/${seg}.mafft.log"
done

# ---------- concatenate by sample name (robust) ----------
echo "==> Concatenating and building single partitioned tree ..."

# Sample order file (use SAMPLES array order)
ORDER="$ALN_DIR/order.txt"; : > "$ORDER"
for sid in "${SAMPLES[@]}"; do echo "$sid" >> "$ORDER"; done

# Build per-segment TSV: "sample<TAB>aligned_sequence" (strip "|SEGMENT" to sample)
declare -a SEG_LEN SEG_TSV
idx=0
for seg in "${SEG_ARR[@]}"; do
  TSV="$ALN_DIR/${seg}.tsv"
  SEG_TSV[$idx]="$TSV"
  awk '
    BEGIN{RS=">"; FS="\n"}
    NR>1{
      header=$1; sub(/\r$/,"",header);
      name=header; sub(/\|.*/,"",name);
      seq="";
      for(i=2;i<=NF;i++){ s=$i; gsub(/[[:space:]]/,"",s); seq=seq s }
      if(length(seq)>0){ print name "\t" seq }
    }' "$ALN_DIR/${seg}.aln.fa" \
    | sort > "$TSV"
  # segment aligned length (assume constant across samples)
  SEG_LEN[$idx]=$(awk -F'\t' 'NR==1{print length($2); exit}' "$TSV")
  ((idx++))
done

# Build concatenated alignment by matching sample names across TSVs
CONCAT="$ALN_DIR/concat.aln.fa"
awk -v OFS="\t" \
    -v t1="${SEG_TSV[0]}" -v t2="${SEG_TSV[1]}" -v t3="${SEG_TSV[2]}" \
    -v L1="${SEG_LEN[0]}" -v L2="${SEG_LEN[1]}" -v L3="${SEG_LEN[2]}" \
    '
    function Ns(n,   s){ s=""; for(i=0;i<n;i++) s=s "N"; return s }
    BEGIN{
      while( (getline < t1) > 0 ){ split($0,a,"\t"); A[a[1]]=a[2] }
      while( (getline < t2) > 0 ){ split($0,a,"\t"); B[a[1]]=a[2] }
      while( (getline < t3) > 0 ){ split($0,a,"\t"); C[a[1]]=a[2] }
    }
    FNR==NR { order[++n]=$1; next } # will be filled by second awk invocation
    END{
      for(i=1;i<=n;i++){
        name=order[i]
        s1=(name in A)?A[name]:Ns(L1)
        s2=(name in B)?B[name]:Ns(L2)
        s3=(name in C)?C[name]:Ns(L3)
        printf(">%s\n%s%s%s\n", name, s1, s2, s3)
      }
    }' "$ORDER" "$ORDER" > "$CONCAT"

# partitions
start=1; PART="$ALN_DIR/concat.partitions"; : > "$PART"
for i in 0 1 2; do
  len="${SEG_LEN[$i]}"; seg="${SEG_ARR[$i]}"
  end=$((start + len - 1))
  echo "DNA, $seg = ${start}-${end}" >> "$PART"
  start=$((end + 1))
done

# diagnostics
nseq=$(alignment_nseq "$CONCAT")
uniq=$(alignment_unique "$CONCAT")
vs=$(alignment_varsites "$CONCAT")
echo "DIAG: concat.aln.fa  nseq=$nseq  unique=$uniq  varSites=$vs"

# ---------- tree from concatenated alignment ----------
TREE_DIR="$OUTDIR/trees"; mkdir -p "$TREE_DIR"
CONCAT_LOG="$OUTDIR/logs/concat.iqtree.log"

if [[ "${uniq:-0}" -gt 1 && "${vs:-0}" -gt 0 ]]; then
  if [[ -n "$IQPROG" ]]; then
    if "$IQPROG" -s "$CONCAT" $PART_FLAG "$PART" -m "$MODEL_FLAG" -B "$BOOT" -alrt "$BOOT" $THREAD_FLAG -pre "$TREE_DIR/concat" > "$CONCAT_LOG" 2>&1; then
      :
    else
      echo "WARN: IQ-TREE failed; trying safer GTR+G, 1 thread" | tee -a "$CONCAT_LOG"
      if "$IQPROG" -s "$CONCAT" $PART_FLAG "$PART" -m GTR+G -B "$BOOT" -alrt "$BOOT" $SAFE_THREAD_FLAG -pre "$TREE_DIR/concat.gtrg" >> "$CONCAT_LOG" 2>&1; then
        mv "$TREE_DIR/concat.gtrg.treefile" "$TREE_DIR/concat.treefile" 2>/dev/null || true
      elif [[ -n "$FASTTREE" ]]; then
        echo "WARN: Falling back to FastTree with bootstrap" | tee -a "$CONCAT_LOG"
        "$FASTTREE" -nt -gtr -gamma -boot "$BOOT" < "$CONCAT" > "$TREE_DIR/concat.treefile" 2>>"$CONCAT_LOG"
      else
        echo "WARN: All ML attempts failed; writing star tree." | tee -a "$CONCAT_LOG"
        write_star_tree "$CONCAT" "$TREE_DIR/concat.treefile"
      fi
    fi
  else
    echo "INFO: Using FastTree with bootstrap" | tee -a "$CONCAT_LOG"
    "$FASTTREE" -nt -gtr -gamma -boot "$BOOT" < "$CONCAT" > "$TREE_DIR/concat.treefile" 2>>"$CONCAT_LOG"
  fi
else
  echo "INFO: Concatenated alignment invariant; writing star tree." | tee -a "$CONCAT_LOG"
  write_star_tree "$CONCAT" "$TREE_DIR/concat.treefile"
fi

# ---------- report ----------
echo
echo "==> Done."
echo "Segments used: ${SEG_ARR[*]}"
echo "Single combined tree:"
echo "  $TREE_DIR/concat.treefile"
echo "Supporting files:"
echo "  Alignment : $ALN_DIR/concat.aln.fa"
echo "  Partitions: $ALN_DIR/concat.partitions"
echo "Logs       : $CONCAT_LOG"
echo
