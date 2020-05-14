#!/bin/bash
#SBATCH --workdir=/home/hugot/projects/20190417_seedgerm_ColxNo_f2
#SBATCH -o scripts/shell/04_snp_naive_freebayes.stdout
#SBATCH -n 1
#SBATCH --mem=30G
#SBATCH --job-name=SNP_calling

# sbatch ~/projects/20190417_seedgerm_ColxNo_f2/scripts/shell/04_snp_naive_freebayes.sh

### Setup

# make output directory
mkdir -p analysis/snp_call/

# reference genome
REF="/home/hugot/reference/arabidopsis/tair10/ensembl_release-37/genome/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa"


### Prepare Freebayes files

# run freebayes naive call
freebayes \
  -f $REF \
  --bam ./mapping/ColxNo_bulk_E1/ColxNo_bulk_E1.bwa.sort.reorder.markdup.rlgn.bam \
  --bam ./mapping/ColxNo_bulk_L1/ColxNo_bulk_L1.bwa.sort.reorder.markdup.rlgn.bam \
  --bam ./mapping/ColxNo_bulk_L2/ColxNo_bulk_L2.bwa.sort.reorder.markdup.rlgn.bam \
  --pooled-continuous \
  --no-complex --no-indels --haplotype-length 5 \
  --min-mapping-quality 20 --min-base-quality 20 \
  --min-alternate-fraction 0.01 --min-alternate-count 2 \
  --min-coverage 10 --max-coverage 400 |\
  bgzip --index --index-name ./analysis/snp_call/ColxNo_all.vcf.gz.gzi -c \
  > ./analysis/snp_call/ColxNo_all_naive.vcf.gz

# Index file
gatk-launch --java-options "-Xms25g -Xmx25g" IndexFeatureFile \
  -F ./analysis/snp_call/ColxNo_all_naive.vcf.gz

# Extract variants to table
gatk-launch --java-options "-Xms25g -Xmx25g" VariantsToTable \
  -V ./analysis/snp_call/ColxNo_all_naive.vcf.gz \
  -F CHROM -F POS -F REF -F ALT -F QUAL -GF AD -GF DP -GF GQ -GF GT \
  -O ./analysis/snp_call/ColxNo_all_naive.tsv



