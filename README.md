Nextflow pipeline for mapping Oxford Nanopore sequencing reads to a reference genome using minimap2. The workflow sorts BAM files, calculates coverage statistics per sample using samtools, and organizes outputs in a results/ directory for easy downstream analysis. Designed for rapid, parallel mapping of multiple samples from a directory of FASTQ files.
output from this are saved in the folder called results
To perform the nextflow mapping execute it like this
nextflow run read_mapping.nf --reads /home/malarm/stx_minimap2_alignment/reads --reference /home/malarm/stx_minimap2_alignment/reference.fasta

