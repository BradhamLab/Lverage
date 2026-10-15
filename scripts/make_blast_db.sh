#!/usr/bin/env bash
# Build the local protein BLAST database that Lverage V2 searches for orthologs.
#
# Give it the NCBI taxonomy IDs of the species you want orthologs from. It
# downloads each species' RefSeq reference proteome from NCBI, combines them,
# builds a BLAST database with -parse_seqids (required: Lverage fetches each
# hit's complete protein with "blastdbcmd -entry"), and verifies retrieval.
#
# Usage:
#   scripts/make_blast_db.sh -t 9606,10090,7227 -o databases/orthologs
#
# Other sources:
#   -i FILE     build from a protein FASTA you already have, instead of
#               downloading. Headers must carry the species in brackets,
#               e.g. ">NP_000001.1 Homeobox protein [Homo sapiens]".
#   -n DB       extract -t TAXIDS from an existing taxonomy-enabled BLAST
#               database (such as NCBI nr) instead of downloading.
#
# Options:
#   -t TAXIDS   comma-separated NCBI taxonomy IDs (find them at
#               https://www.ncbi.nlm.nih.gov/taxonomy)
#   -o PREFIX   output database prefix, e.g. databases/orthologs (required)
#   -h          show this help
#
# Requires BLAST+ (makeblastdb, blastdbcmd) on PATH. Downloading also requires
# the NCBI datasets CLI and unzip. With conda:
#   conda install -c conda-forge -c bioconda blast ncbi-datasets-cli unzip

set -euo pipefail

usage() { sed -n '2,/^set -euo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; }
die() { echo "error: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null || die "$1 not found on PATH. $2"; }

input=""; prefix=""; source_db=""; taxids=""
while getopts "i:o:n:t:h" opt; do
    case "$opt" in
        i) input="$OPTARG" ;;
        o) prefix="$OPTARG" ;;
        n) source_db="$OPTARG" ;;
        t) taxids="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage >&2; exit 2 ;;
    esac
done

[[ -n "$prefix" ]] || die "-o PREFIX is required"
need makeblastdb "Install BLAST+ first."
need blastdbcmd "Install BLAST+ first."
mkdir -p "$(dirname "$prefix")"

if [[ -n "$input" && ( -n "$source_db" || -n "$taxids" ) ]]; then
    die "use -i FILE on its own, or -t TAXIDS (optionally with -n DB)"
fi

if [[ -n "$source_db" ]]; then
    [[ -n "$taxids" ]] || die "-t TAXIDS is required with -n"
    input="${prefix}.fasta"
    echo "Extracting taxids $taxids from $source_db into $input"
    # -target_only keeps only deflines belonging to the requested taxa, so the
    # species in each header is one Lverage's species map will recognise.
    blastdbcmd -db "$source_db" -taxids "$taxids" -target_only -outfmt '%f' -out "$input"

elif [[ -n "$taxids" ]]; then
    need datasets "Install the NCBI datasets CLI (conda: bioconda::ncbi-datasets-cli)."
    need unzip "Install unzip."
    input="${prefix}.fasta"
    : > "$input"
    workdir=$(mktemp -d "${TMPDIR:-/tmp}/make_blast_db.XXXXXX")
    trap 'rm -rf "$workdir"' EXIT
    IFS=',' read -r -a taxid_list <<< "$taxids"
    for taxid in "${taxid_list[@]}"; do
        taxid="${taxid//[[:space:]]/}"
        [[ "$taxid" =~ ^[0-9]+$ ]] || die "taxid is not a number: '$taxid'"
        echo "Downloading RefSeq reference proteome for taxid $taxid"
        datasets download genome taxon "$taxid" --reference --include protein --no-progressbar \
            --filename "$workdir/$taxid.zip" \
            || die "download failed for taxid $taxid (no RefSeq reference genome, or no network?)"
        unzip -q -o "$workdir/$taxid.zip" -d "$workdir/$taxid"
        faa=$(find "$workdir/$taxid" -name 'protein.faa' | head -n 1)
        [[ -n "$faa" ]] || die "no protein.faa in the download for taxid $taxid; its reference genome may lack annotation"
        count=$(grep -c '^>' "$faa")
        species=$(grep -m1 '^>' "$faa" | sed -n 's/.*\[\(.*\)\].*/\1/p')
        echo "  $count proteins, species in headers: ${species:-NONE}"
        [[ -n "$species" ]] || die "headers for taxid $taxid have no bracketed species name; Lverage cannot use them"
        # RefSeq reuses WP_ accessions for identical proteins across related
        # bacteria, and makeblastdb -parse_seqids refuses duplicate ids. When a
        # record's id was already added by an earlier species, suffix it with
        # this taxid; the description and species are kept as they are.
        { grep '^>' "$input" || true; } | sed 's/^>//; s/[[:space:]].*//' > "$workdir/seen_ids"
        awk -v taxid="$taxid" -v seen="$workdir/seen_ids" '
            BEGIN { while ((getline id < seen) > 0) have[id] = 1 }
            /^>/ { id = substr($1, 2)
                   if (id in have) { $1 = ">" id "_" taxid; renamed++ }
                   print; next }
            { print }
            END { if (renamed) printf("  %d ids were already present from another species and were suffixed with _%s\n", renamed, taxid) > "/dev/stderr" }
        ' "$faa" >> "$input"
    done
    echo "Combined proteome written to $input"
fi

[[ -n "$input" ]] || die "give -t TAXIDS to download, or -i FILE to build from your own FASTA"
[[ -s "$input" ]] || die "input FASTA is missing or empty: $input"

dup_ids=$(grep '^>' "$input" | sed 's/^>//; s/[[:space:]].*//' | sort | uniq -d | sed -n '1,5p')
if [[ -n "$dup_ids" ]]; then
    die "duplicate record identifiers in $input (makeblastdb -parse_seqids requires unique ids); first few: $(echo "$dup_ids" | tr '\n' ' ')"
fi

echo "Building $prefix"
makeblastdb -in "$input" -dbtype prot -parse_seqids -max_file_sz 3GB -out "$prefix" >/dev/null

# Verify that accession retrieval works, using the first record in the input.
first_id=$(grep -m1 '^>' "$input" | sed 's/^>//; s/[[:space:]].*//')
[[ -n "$first_id" ]] || die "could not read a record identifier from $input"
blastdbcmd -db "$prefix" -entry "$first_id" -outfmt '%a' >/dev/null \
    || die "built $prefix but could not retrieve $first_id from it; Lverage will not be able to use this database"

echo
blastdbcmd -db "$prefix" -info | sed -n '1,2p'
echo
echo "Species found in headers (use these exact names as species_map keys):"
grep '^>' "$input" | sed -n 's/.*\[\([^]]*\)\].*/\1/p' | sort | uniq -c | sort -rn
echo
echo "Database ready. Pass database_path=\"$prefix\" to LocalBlastSearcher."
