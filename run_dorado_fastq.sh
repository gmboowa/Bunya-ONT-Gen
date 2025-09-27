#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------
# Usage:
#   ./run_dorado_fastq.sh -i <input_dir> -o <output_dir> -s <sample_sheet.csv> -k <kit_name> [-- extra dorado args...]
# ---------------------------------------

DORADO_BIN="${DORADO_BIN:-dorado}"

IN_DIR=""
OUT_DIR=""
SAMPLE_SHEET=""
KIT_NAME=""

print_usage() {
  cat >&2 <<'EOF'
Usage: run_dorado_fastq.sh -i <input_dir> -o <output_dir> -s <sample_sheet.csv> -k <kit_name> [-- extra dorado args...]
Required:
  -i  Input directory containing pooled *.fastq.gz / *.fq.gz / *.fastq / *.fq
  -o  Output directory root (demuxed FASTQs will go to <output_dir>/demux/)
  -s  Sample sheet CSV (per Dorado format)
  -k  Kit name (e.g., TWIST-96A-UDI)
Notes:
  - Script forces FASTQ output by adding '--emit-fastq' to 'dorado demux'.
  - Any args after -- are appended to 'dorado demux' (you may override defaults).
Env:
  DORADO_BIN  Path to dorado binary (default: 'dorado' in PATH)
EOF
}

# ---- Parse flags ----
while getopts ":i:o:s:k:h" opt; do
  case "$opt" in
    i) IN_DIR="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    s) SAMPLE_SHEET="$OPTARG" ;;
    k) KIT_NAME="$OPTARG" ;;
    h) print_usage; exit 0 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; print_usage; exit 2 ;;
    :)  echo "Missing argument for -$OPTARG" >&2; print_usage; exit 2 ;;
  esac
done
shift $((OPTIND-1))

# Capture extra args after --
EXTRA_ARGS=()   # init to avoid "unbound variable" with set -u
if [[ "${1:-}" == "--" ]]; then
  shift
  # shellcheck disable=SC2206
  EXTRA_ARGS=($@)
fi

# ---- Validations ----
[[ -n "$IN_DIR" && -d "$IN_DIR" ]] || { echo "ERROR: -i <input_dir> missing or invalid." >&2; exit 1; }
[[ -n "$OUT_DIR" ]] || { echo "ERROR: -o <output_dir> is required." >&2; exit 1; }
[[ -n "$SAMPLE_SHEET" && -f "$SAMPLE_SHEET" ]] || { echo "ERROR: -s <sample_sheet.csv> missing or invalid." >&2; exit 1; }
[[ -n "$KIT_NAME" ]] || { echo "ERROR: -k <kit_name> is required." >&2; exit 1; }

if ! command -v "$DORADO_BIN" >/dev/null 2>&1; then
  echo "ERROR: dorado not found. Set DORADO_BIN=/path/to/dorado or add to PATH." >&2
  exit 1
fi

# Root output folder where all FASTQs end up (flat)
demux_root="$OUT_DIR/demux"
mkdir -p "$demux_root"
echo "[info] All outputs will go to: $demux_root"

# Iterate pooled FASTQs
found_any=false
while IFS= read -r -d '' fq; do
  found_any=true
  in_base="$(basename "$fq")"

  # Input stem (strip .gz and .fastq/.fq)
  stem="$in_base"
  stem="${stem%.gz}"
  stem="${stem%.fastq}"
  stem="${stem%.fq}"

  # Per-input staging dir; dorado writes subfolders here
  stage_dir="$demux_root/.stage_${stem}_$$"
  mkdir -p "$stage_dir"

  echo "[info] Demuxing: $in_base -> $stage_dir"
  set -x
  "$DORADO_BIN" demux \
    --sample-sheet "$SAMPLE_SHEET" \
    --kit-name "$KIT_NAME" \
    --emit-fastq \
    --output-dir "$stage_dir" \
    ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"} \
    "$fq"
  { set +x; } 2>/dev/null

  # Recursively collect FASTQs from dorado's subfolders and move into flat demux_root
  moved_any=false
  while IFS= read -r -d '' f; do
    moved_any=true
    base="$(basename "$f")"
    # Ensure gzipped output in the flat folder, prefix with input stem to avoid collisions
    if [[ "$base" != *.gz ]]; then
      gzip -c "$f" > "$demux_root/${stem}__${base}.gz"
    else
      mv -f "$f" "$demux_root/${stem}__${base}"
    fi
  done < <(find "$stage_dir" -type f \( -name '*.fastq' -o -name '*.fq' -o -name '*.fastq.gz' -o -name '*.fq.gz' \) -print0)

  if ! $moved_any; then
    echo "[warn] No FASTQ files found under: $stage_dir"
    echo "[warn] Dorado may have only written BAM (classification)."
    echo "[warn] Contents:"
    find "$stage_dir" -maxdepth 3 -print | sed 's/^/[warn]   /'
  fi

  # Clean up staging dir
  rm -rf "$stage_dir" 2>/dev/null || true

done < <(find "$IN_DIR" -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' -o -name '*.fastq' -o -name '*.fq' \) -print0)

if ! $found_any; then
  echo "ERROR: No FASTQ files found under: $IN_DIR" >&2
  exit 1
fi

echo "[info] Done. Demultiplexed files are in: $demux_root"
