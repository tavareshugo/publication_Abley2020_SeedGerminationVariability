#!/bin/bash
#SBATCH --job-name=seed_filter
#SBATCH -o /home/hugot/projects/20190417_seedgerm_ColxNo_f2/scripts/shell/01_quality_filtering.stdout
#SBATCH -n 16
#SBATCH --mem-per-cpu=1G

## sbatch $HOME/projects/20190417_seedgerm_ColxNo_f2/scripts/shell/01_quality_filtering.sh

##########
# Compile sample information and quality filter
##########

wd="$HOME/projects/20190417_seedgerm_ColxNo_f2/"
fastq_dir="$HOME/sequence_data/DNAseq/2019_poolseq_seedvariability/X201SC19030178-Z01-F001/raw_data"
r1_suffix="_1.fq.gz"
r2_suffix="_2.fq.gz"

#
# Prepare output directory
#
cd ${wd}
mkdir -p reads/


#
# Make a list of file name prefixes
#
find ${fastq_dir} -type f -name "*${r1_suffix}" | \
sed "s/${r1_suffix}//" > ./reads/raw_read_file_prefix.txt


#
# Extract sample information
#
## Create header of output file
printf "file_path,file_prefix,sample_name,barcode_seq,run_id\n" > ./reads/read_info.csv

while read file
do
  # get path to read files
  read_path=$(dirname $file)

  # get file prefix basename
  read_basename=$(basename $file)

  # sample name in this case is the directory name
  sample_name=$(dirname $file | sed 's/.*\///g')

  # Get barcode from read header
  #barcode_seq=$(zcat ${file}${r1_suffix} | head -n 1 | cut -d ' ' -f2 | sed 's/.*://g')
  barcode_seq=$(zcat ${file}${r1_suffix} | sed -n '1~4p' | head -n 500 | sed 's/.*://' | sort | uniq -c | sort -nr | head -n 1 | sed 's/.* //')

  # Get run id (machine ID + Run ID) from read header
  run_id=$(zcat ${file}${r1_suffix} | head -n 1 | cut -d ' ' -f1 | awk -F ":" '{print $1 ":" $2 ":" $3}')

  # Add result to output file
  printf "$read_path,$read_basename,$sample_name,$barcode_seq,$run_id\n" >> ./reads/read_info.csv

done < ./reads/raw_read_file_prefix.txt
wait


#
# Filtering and QC
#
# script available from https://github.com/tavareshugo/BioPipelines/blob/master/qc_fastq.sh
# Remove reads with any N
# trim based on quality 20
# keep reads with at least 50bp length
# cut the last nucleotide of the read (cycle 151, which usually has poor quality)
$HOME/code/bioPipelines/qc_fastq.sh -c $SLURM_NPROCS \
-o "${wd}/reads/" \
-p "${wd}/reads/raw_read_file_prefix.txt" \
-1 "${r1_suffix}" -2 "${r2_suffix}" \
-f "--trim-n --quality-cutoff 20 --minimum-length 50 --max-n 0 --cut -1"
