#!/bin/bash

###############
# Aligning reads to genome
###############

# This submits a load of jobs to the scheduler. Run as:
# bash $HOME/projects/20190417_seedgerm_ColxNo_f2/scripts/shell/02_mapping.sh

#
# Adjust these paths as necessary ####
#
#working directory
wd="$HOME/projects/20190417_seedgerm_ColxNo_f2"

#reference genome
ref="$HOME/reference/arabidopsis/tair10/ensembl_release-37/genome/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa"


#
# Set things up
#
# Change to working directory
cd "${wd}"

# Make directory for output
mkdir -p mapping/
out_dir=${wd}/mapping

# Path to filtered reads
read_path="${wd}/reads/filtered_reads/"

# Get fastq file prefix
reads=$(find ${read_path} -type f -name *_1.fq.gz | sed 's/_1.fq.gz//')

# Suffix
r1_suffix="_1.fq.gz"
r2_suffix="_2.fq.gz"


for sample in $reads
do
  read1=${sample}${r1_suffix}
  read2=${sample}${r2_suffix}

  # This big pipe takes the first few reads of the fastq file and gets the most
  ## commonly represented barcode annotated in the file
  barcode=$(zcat ${read1} | head -n 1000 | grep @ | cut -d " " -f 2 | sed 's/.*://' | sort | uniq -c | sort -n | sed 's/.* //' | tail -n 1)

  # Sample name
  sample_name=$(basename ${sample} | sed 's/_.*//')
  sample_name="ColxNo_bulk_${sample_name}" # add informative prefix

  # make output directory
  mkdir -p ${out_dir}/${sample_name}

  # Launch job
  echo \
  '#!/bin/bash
  /home/hugot/code/bioPipelines/mapping_pipeline.sh \
  -m bwa \
  -o ${out_dir}/$sample_name/ \
  -p $sample_name \
  -1 $read1 \
  -2 $read2 \
  -i $sample_name \
  -l $sample_name \
  -b $barcode \
  -s $sample_name \
  -r $ref \
  -c 16 \
  -d yes \
  -k yes' | \
  sbatch --job-name=${sample_name}_map \
  -o $out_dir/${sample_name}/${sample_name}.stdout \
  --export=ALL \
  -n 16 --workdir=${wd}
done


# # Run multiqc once mapping is over
# sbatch --job-name=multiqc_seed \
# -o ${out_dir}/multiqc.stdout \
# --wrap "multiqc --filename mapping_reports.html --outdir ${out_dir} ${out_dir}/*/stats"
