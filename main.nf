#!/usr/bin/env nextflow
params.name             = "Rinn Lab RNA-seq test"
params.email            = "michael.smallegan@colorado.edu"
params.reads            = "${baseDir}/input/fastq/*{*_R1,*_R2}.fastq.gz"
params.genome           = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/GRCm38.primary_assembly.genome.fa.gz"
params.annotation       = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gtf.gz"



log.info "Rinn Lab RNA-seq Pipeline"
log.info "====================================="
log.info "name         : ${params.name}"
log.info "email        : ${params.email}"
log.info "reads        : ${params.reads}"
log.info "genome       : ${params.genome}"
log.info "annotation   : ${params.annotation}"
log.info "\n"



reads = Channel.fromFilePairs(params.reads, size: -1)
  .ifEmpty { error "Can't find any reads matching: ${params.reads}" }
  .into {
    reads_for_fastqc;
    reads_for_mapping
  }


process retrieve_annotation {

  publishDir 'results/annotation'

  output:
  file "*.gtf" into annotation

  script:
  """
  wget ${params.annotation}
  gunzip *.gz
  """
}



annotation.into {
  annotation_for_index;
  annotation_for_count;
  annotation_for_transcriptome
}



process retrieve_genome {

  publishDir 'results/genome'

  output:
  file('*.fa') into genome

  script:
  """
  wget ${params.genome}
  gunzip *.gz
  """
}



genome.into {
  genome_for_transcriptome;
  genome_for_index
}



process make_transcriptome {

  publishDir 'results/genome'

  input:
  file annotation from annotation_for_transcriptome
  file genome from genome_for_transcriptome

  output:
  file "transcriptome.fa" into transcriptome

  script:
  """
  gffread -w transcriptome.fa -g ${genome} ${annotation}
  """
}



process fastqc {

  publishDir 'results/fastqc'

  input:
  set sample_id, file(fastqz) from reads_for_fastqc

  output:
  file "*.zip" into fastqc

  script:
  """
  fastqc --threads 4 -f fastq -q ${fastqz}
  """
}



process index {

  memory '60 GB'
  clusterOptions '--nodes=1 --ntasks=64'
  publishDir 'results/star_index'

  input:
  file genome_for_index
  file annotation_for_index

  output:
  file "index" into index

  script:
  """
  mkdir index
  STAR  --runThreadN 64 \
  --runMode genomeGenerate \
  --genomeDir ./index \
  --genomeFastaFiles ${genome_for_index} \
  --sjdbGTFfile ${annotation_for_index}
  """
}



process map {

  clusterOptions '--nodes=1 --ntasks=4 --mem=40gb'
  publishDir 'results/bam'

  input:
  set sample_id, file(reads), file(index) from reads_for_mapping.combine(index)

  output:
  set sample_id, file("*Aligned.out.bam") into mapped_genome
  set sample_id, file("*toTranscriptome.out.bam") into mapped_transcriptome
  file '*' into star

  script:
  """
  STAR  --runThreadN 4 \
  --genomeDir ${index} \
  --readFilesIn ${reads.findAll{ it =~ /\_R1\./ }.join(',')} \
                ${reads.findAll{ it =~ /\_R2\./ }.join(',')} \
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted \
  --outSAMunmapped Within \
  --outSAMattributes NH HI NM MD AS \
  --outReadsUnmapped Fastx \
  --quantMode TranscriptomeSAM \
  --outFileNamePrefix ${sample_id}_
  """
}



mapped_genome.into {
  mapped_for_count;
  mapped_for_igv;
  mapped_for_markduplicates
}



process count {

  publishDir 'results/feature_counts'

  input:
  file annotation from annotation_for_count
  set sample_id, file(bam) from mapped_for_count

  output:
  file '*.fCounts'
  file '*.fCounts*' into counts

  script:
  """
  featureCounts  -C \
    -p \
    -T 4 \
    -g gene_id \
    -a ${annotation} \
    -o ${sample_id}.fCounts \
    ${bam}
  """
}



process salmon {

  publishDir 'results/salmon'
  clusterOptions '--ntasks=8'

  input:
  file transcript_fasta from transcriptome
  set sample_id, file(bam) from mapped_transcriptome

  output:
  file("*") into salmon

  script:
  """
  salmon quant -l A \
    -p 8 \
    -t ${transcript_fasta} \
    -o ${sample_id} \
    -a ${bam} \
    --numBootstraps 30
  """
}



process sort_bam {

  publishDir 'results/igv'
  clusterOptions '--mem=32gb'

  input:
  set sample_id, file(bam_file) from mapped_for_igv

  output:
  set sample_id, file("*.bam"), file('*.bai') into sorted_bam

  script:
  """
  samtools sort --threads 4 \
    -m 4G \
    -o ${sample_id}.bam \
    ${bam_file}
  samtools index ${sample_id}.bam
  """
}



process compile_fastqc {

  input:
  file fastqc from fastqc.collect()

  output:
  file "fastqc" into fastqc_compiled

  script:
  """
  mkdir fastqc
  mv ${fastqc} fastqc/.
  """
}



process compile_star {

  input:
  file star from star.collect()

  output:
  file "star" into star_compiled

  script:
  """
  mkdir star
  mv ${star} star/.
  """
}



process compile_salmon {

  input:
  file salmon from salmon.collect()

  output:
  file "salmon" into salmon_compiled

  script:
  """
  mkdir salmon
  mv ${salmon} salmon/.
  """
}



process compile_counts {

  input:
  file counts from counts.collect()

  output:
  file "counts" into counts_compiled

  script:
  """
  mkdir counts
  mv ${counts} counts/.
  """
}



process multiqc {

  publishDir "results/multiqc"
  clusterOptions '--ntasks=1'

  input:
  file fastqc from fastqc_compiled
  file star from star_compiled
  file salmon from salmon_compiled
  file counts from counts_compiled

  output:
  set file('*_multiqc_report.html'), file('*_data/*')

  script:
  """
  export LC_ALL=en_US.UTF-8
  export LANG=en_US.UTF-8

  multiqc ${fastqc} \
    ${star} \
    ${salmon} \
    ${counts} \
    --title '${params.name}' \
    --cl_config "extra_fn_clean_exts: [ '_1', '_2' ]"
  """
}



// process differential_expression {
//
//   publishDir 'reports', mode: 'copy'
//
//   input:
//   file sample_file from salmon_out.collect()
//
//   script:
//   """
//   Rscript -e 'rmarkdown::render("${baseDir}/analysis/01_quality_control.Rmd")'
//   Rscript -e 'rmarkdown::render("${baseDir}/analysis/02_differential_expression.Rmd")'
//   Rscript -e 'rmarkdown::render("${baseDir}/analysis/03_functional_analysis.Rmd")'
//   """
// }
