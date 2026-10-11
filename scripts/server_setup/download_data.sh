#!/usr/bin/env bash
# Login-node network preparation, with atomic downloads and provenance.
set -euo pipefail
repo=${1:?Repository source required}
source "$(dirname "${BASH_SOURCE[0]}")/environment.sh"
module load BCFtools/1.22-GCC-14.2.0
data="$FMB_CAMPAIGN_ROOT/data"
mkdir -p "$data/legacy-wide-vcf" "$data/native-eur-vcf" "$data/genetic_maps"
cd "$data"
base=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
curl --fail --location --retry 3 "$base/integrated_call_samples_v3.20130502.ALL.panel" -o samples.panel.part
mv samples.panel.part samples.panel
awk 'NR>1 && $3=="EUR" {print $1}' samples.panel > eur-samples.txt
test -s eur-samples.txt
for chrom in $(seq 1 22); do
  file="$data/genetic_maps/genetic_map_GRCh37_chr${chrom}.txt.gz"
  if [[ ! -s "$file" ]]; then
    curl --fail --location --retry 3 "https://raw.githubusercontent.com/adimitromanolakis/geneticMap-GRCh37/master/genetic_map_GRCh37_chr${chrom}.txt.gz" -o "$file.part"
    gzip -t "$file.part"
    mv "$file.part" "$file"
  fi
done
printf 'track\tregion\tchrom\tstart\tend\tpopulation\tvariants\tsource\n' > download-regions.tsv
download_region() {
  local track=$1 rid=$2 chrom=$3 start=$4 end=$5 population=$6 directory=$7
  local max_end=${8:-$end} original_start=$start width=$((end-start)) count=0
  local url="$base/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  local output="$directory/$rid.vcf.gz"
  if [[ -s "$output.window.tsv" ]]; then
    read -r start end < "$output.window.tsv"
  fi
  if [[ -s "$output" && -s "$output.tbi" ]]; then count=$(bcftools index -n "$output"); fi
  while [[ "$count" -le 0 ]]; do
    if [[ "$population" == EUR ]]; then
      tabix -h "$url" "$chrom:$start-$end" | bcftools view -S "$data/eur-samples.txt" -m2 -M2 -v snps -Oz -o "$output.part"
    else
      tabix -h "$url" "$chrom:$start-$end" | bcftools view -m2 -M2 -v snps -Oz -o "$output.part"
    fi
    gzip -t "$output.part"
    count=$(bcftools view -H "$output.part" | wc -l)
    if [[ "$count" -le 0 ]]; then
      printf 'Empty source window rejected: %s chr%s:%s-%s\n' "$rid" "$chrom" "$start" "$end" >&2
      rm -f "$output.part" "$output.part.tbi"
      [[ "$track" == native_subset ]] || exit 1
      start=$((start+1000000)); end=$((start+width))
      if [[ "$end" -ge "$max_end" ]]; then
        printf 'No accessible native window before neighbouring region: %s\n' "$rid" >&2
        exit 1
      fi
      continue
    fi
    tabix -p vcf "$output.part"
    mv "$output.part" "$output"
    mv "$output.part.tbi" "$output.tbi"
  done
  printf '%s\t%s\n' "$start" "$end" > "$output.window.tsv"
  test "$count" -gt 0
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$track" "$rid" "$chrom" "$start" "$end" "$population" "$count" "$url" >> download-regions.tsv
  if [[ "$track" == native_subset ]]; then
    printf '%s\t%s\t%s\t%s\t%s\n' "$rid" "$chrom" "$start" "$end" "$original_start" >> native-windows.tsv
  fi
}
# New wide windows preserve declared SNP-size capacity. They are explicitly a
# new data realization, not the inaccessible historical VCF/prediction archive.
while IFS=, read -r rid chrom start end notes; do
  [[ "$rid" == region_id ]] && continue
  download_region legacy_reconstruction "$rid" "$chrom" "$start" "$((start+2000000))" ALL "$data/legacy-wide-vcf"
done < "$repo/inst/extdata/regions.csv"
python3 - "$repo/inst/extdata/gwfm_regions.csv" > native-candidates.tsv <<'PY'
import csv, sys
with open(sys.argv[1]) as src:
    rows = list(csv.DictReader(src))
for row in rows:
    neighbours = [int(x['start']) for x in rows if x['chrom'] == row['chrom'] and int(x['start']) > int(row['start'])]
    boundary = min(neighbours) if neighbours else int(row['start']) + 10000000
    print(row['region_id'], row['chrom'], row['start'], row['end'], boundary, sep='\t')
PY
printf 'region_id\tchrom\tstart\tend\toriginal_start\n' > native-windows.tsv
while IFS=$'\t' read -r rid chrom start end max_end; do
  download_region native_subset "$rid" "$chrom" "$start" "$end" EUR "$data/native-eur-vcf" "$max_end"
done < native-candidates.tsv
python3 - "$data" <<'PY'
import csv, pathlib, sys
root = pathlib.Path(sys.argv[1])
with (root/'native-windows.tsv').open() as src, (root/'native-regions.csv').open('w') as dst:
    writer = csv.DictWriter(dst, fieldnames=['region_id','chrom','start','end','notes'])
    writer.writeheader()
    for row in csv.DictReader(src, delimiter='\t'):
        notes = 'EUR source window; original start '+row.pop('original_start')
        writer.writerow(dict(row, notes=notes))
PY
find "$data" -type f ! -name SHA256SUMS ! -name inputs-manifest.csv ! -name DATA_COMPLETE ! -name '*.part' -print0 | sort -z | xargs -0 sha256sum > "$data/SHA256SUMS"
python3 - "$data" <<'PY'
import csv, pathlib, sys
root = pathlib.Path(sys.argv[1])
with (root / "SHA256SUMS").open() as src, (root / "inputs-manifest.csv").open("w") as dst:
    writer = csv.writer(dst)
    writer.writerow(["path", "sha256"])
    for line in src:
        digest, path = line.rstrip("\n").split("  ", 1)
        writer.writerow([path, digest])
PY
touch "$data/DATA_COMPLETE"
