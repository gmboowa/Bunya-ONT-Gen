#!/usr/bin/env bash
set -euo pipefail

# Requirements: snpEff (v5.x), entrez-direct (esearch/efetch)

DB=""; NAME=""; QUERY=""; ACCS=""
usage() {
  cat <<EOF
Usage:
  $0 --db <id> [--name "Human name"] (--query "<Entrez query>" | --acc "ACC1,ACC2,...")

Examples:
  $0 --db cchfv_IbAr10200 --name "CCHFV IbAr10200" \
     --query '"Crimean-Congo hemorrhagic fever virus"[Organism] AND IbAr10200 AND (segment[Title] OR L[Title] OR M[Title] OR S[Title]) AND complete'

  $0 --db rvfv_ZH501 --name "Rift Valley fever virus ZH501" \
     --acc "NC_014395.1,NC_014396.1,NC_014397.1"
EOF
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --db) DB="$2"; shift 2;;
    --name) NAME="$2"; shift 2;;
    --query) QUERY="$2"; shift 2;;
    --acc) ACCS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done
[[ -z "$DB" ]] && { echo "--db is required"; usage; exit 1; }
[[ -z "$QUERY" && -z "$ACCS" ]] && { echo "Provide --query or --acc"; usage; exit 1; }

# Locate / prepare SnpEff home
: "${CONDA_PREFIX:=}"
SNPEFF_SRC="$(ls -d "$CONDA_PREFIX"/share/snpeff-* 2>/dev/null | head -n1 || true)"
[[ -z "$SNPEFF_SRC" ]] && { echo "Cannot find conda snpEff share dir. Activate env with snpEff."; exit 1; }
SNPEFF_HOME="${SNPEFF_HOME:-$HOME/snpeff}"
[[ -d "$SNPEFF_HOME" ]] || cp -R "$SNPEFF_SRC" "$SNPEFF_HOME"

# Register genome name
NAME="${NAME:-$DB}"
grep -q "^$DB\.genome" "$SNPEFF_HOME/snpEff.config" || echo "$DB.genome : $NAME" >> "$SNPEFF_HOME/snpEff.config"

# Fetch GenBank(s)
mkdir -p "$SNPEFF_HOME/data/$DB"
GB="$SNPEFF_HOME/data/$DB/genes.gb"
: > "$GB"

if [[ -n "$ACCS" ]]; then
  IFS=',' read -r -a ARR <<<"$ACCS"
  for A in "${ARR[@]}"; do
    echo "Fetching $A ..."
    efetch -db nucleotide -format gb -id "$A" >> "$GB"
  done
else
  echo "Searching NCBI with query: $QUERY"
  esearch -db nucleotide -query "$QUERY" | efetch -format gb > "$GB"
fi

echo "Records in genes.gb: $(grep -c '^LOCUS' "$GB")"
[[ -s "$GB" ]] || { echo "Empty genes.gb. Check your query/accessions."; exit 1; }

# Build DB
snpEff -c "$SNPEFF_HOME/snpEff.config" build -genbank -v "$DB"
echo "Built DB '$DB'. Annotate with:"
echo "  snpEff -c $SNPEFF_HOME/snpEff.config ann $DB input.vcf > output.annotated.vcf"


