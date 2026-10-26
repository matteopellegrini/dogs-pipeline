#!/bin/bash
# Sample sheet for the SRA breed top-up (panel v2.1): one row per downloaded run that has a
# DONE marker. fastq-dump writes <RUN>_1/_2.fastq.gz; the pipeline wants *_R1_*/*_R2_*, so the
# files are symlinked under those names in place. Rows run stages 1-8 only (TO_STAGE=8:
# align + GLIMPSE imputation, no report), e.g.
#   bash cluster/make-sra-sheet.sh > sample_sheet.sra.tsv
#   qsub -t 2-N -tc 30 cluster/submit-array.sh sample_sheet.sra.tsv 1     (with TO_STAGE=8 exported)
D=/u/project/pellegrini/gkislik/dogs
printf 'sample_id\tfastq_dir\tage\toutput_name\twork_dir\tpub_dir\tfrom_stage\tsex\tnotes\n'
for done in "$D"/fastq/sra/*/DONE; do
  dir=$(dirname "$done"); run=$(basename "$dir"); label=$(cut -d' ' -f1 "$done")
  [[ -s "$dir/${run}_1.fastq.gz" && -s "$dir/${run}_2.fastq.gz" ]] || { echo "WARN $run: missing pair" >&2; continue; }
  ln -sf "${run}_1.fastq.gz" "$dir/${run}_R1_001.fastq.gz"; ln -sf "${run}_2.fastq.gz" "$dir/${run}_R2_001.fastq.gz"
  printf 'SRA-%s\t%s\t\tSRA-%s\t%s/work_sra/%s/analysis\t%s/results_sra/%s\t1\t\t%s\n' "$run" "$dir" "$run" "$D" "$run" "$D" "$run" "$label"
done
