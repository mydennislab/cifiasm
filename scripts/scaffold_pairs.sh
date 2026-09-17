#!/usr/bin/env bash
# Export YaHS contacts as scaffold-coordinate 4DN pairs.
set -euo pipefail

while (( $# )); do
    if (( $# < 2 )); then
        echo "Missing value for $1" >&2
        exit 2
    fi
    case "$1" in
        --bin) bin=$2 ;;
        --agp) agp=$2 ;;
        --ctg-fai) ctg_fai=$2 ;;
        --scaf-fai) scaf_fai=$2 ;;
        --output) pairs=$2 ;;
        --log) log=$2 ;;
        --mapq) mapq=$2 ;;
        --outdir) outdir=$2 ;;
        --sort-tmp) sort_tmp=$2 ;;
        --sort-mem) sort_mem=$2 ;;
        --sort-threads) sort_threads=$2 ;;
        --threads) threads=$2 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
    shift 2
done

: "${bin:?Missing --bin}" "${agp:?Missing --agp}"
: "${ctg_fai:?Missing --ctg-fai}" "${scaf_fai:?Missing --scaf-fai}"
: "${pairs:?Missing --output}" "${log:?Missing --log}"
: "${mapq:?Missing --mapq}" "${outdir:?Missing --outdir}"
: "${sort_tmp:?Missing --sort-tmp}" "${sort_mem:?Missing --sort-mem}"
: "${sort_threads:?Missing --sort-threads}" "${threads:?Missing --threads}"

mkdir -p "$outdir"
rm -rf "$sort_tmp" "$pairs.part"
trap 'rm -rf "$sort_tmp" "$pairs.part"' EXIT
mkdir -p "$sort_tmp"
: > "$log"
{
    printf '## pairs format v1.0\n#sorted: chr1-chr2-pos1-pos2\n#shape: upper triangle\n'
    awk '{ print "#chromsize: " $1 " " $2 }' "$scaf_fai"
    printf '#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n'
    juicer pre -q "$mapq" "$bin" "$agp" "$ctg_fai" 2>> "$log" \
        | awk -v fai="$scaf_fai" '
            BEGIN {
                OFS = "\t"
                while ((getline line < fai) > 0) {
                    split(line, f, "\t")
                    n++
                    rank[f[1]] = n
                    len[f[1]] = f[2]
                }
            }
            {
                c1 = $2; p1 = $3 + 1
                c2 = $6; p2 = $7 + 1
                if (!(c1 in rank) || !(c2 in rank)) {
                    bad_chr++
                    next
                }
                if (p1 < 1 || p2 < 1 || p1 > len[c1] || p2 > len[c2]) {
                    bad_pos++
                    next
                }
                if (rank[c1] > rank[c2] || (rank[c1] == rank[c2] && p1 > p2)) {
                    t = c1; c1 = c2; c2 = t
                    t = p1; p1 = p2; p2 = t
                }
                kept++
                print rank[c1], rank[c2], p1, p2, c1, c2
            }
            END {
                printf "[pairs] %d records written, %d dropped (sequence not in the scaffold fai), %d dropped (position outside the sequence)\n", kept, bad_chr, bad_pos > "/dev/stderr"
            }' 2>> "$log" \
        | LC_ALL=C sort -k1,1n -k2,2n -k3,3n -k4,4n \
            -S "$sort_mem" -T "$sort_tmp" --parallel="$sort_threads" \
        | awk 'BEGIN { OFS = "\t" } { print ".", $5, $3, $6, $4, ".", "." }'
} | bgzip -@ "$threads" > "$pairs.part"

if grep -q 'using scale factor' "$log"; then
    echo "ERROR: juicer pre scaled the coordinates (see $log); the pairs would not be in bp" >&2
    exit 1
fi
mv "$pairs.part" "$pairs"
