# Running deNOVOEnrich on Test samples 
This directory contains test files to validate and demonstrate the usage of the deNOVOEnrich pipeline.

## Note on raw sequencing input files
The paired-end test FASTQ files (Test_1.fq.bz2 & Test_2.fq.bz2) are compressed with bzip2. Before running the pipeline, these files should be extracted, as the pipeline **does not accept .bz2 files** directly.

## Reference genome
The test run uses the Arabidopsis thaliana **Col-CEN reference** assembly.
You can download the FASTA file from TAIR: https://www.arabidopsis.org/download/list?dir=Genes%2FCol-CEN_genome_assembly_release

# Example command to run the test files

<pre> bash deNOVOEnrich.sh --genome ColCEN.fasta --ref_TE . --rawRead1 Test_1.fq --rawRead2 Test_2.fq --Sample Test --adapter /Full/path/to/TruSeq3-PE.fa --outDir ./denovo_Test --TEfam AtCOPIA93 --CORES 10 --monomeric_repeat Col-CEN_monomeric_repeat.bed </pre>

_Please ensure you update all file paths to match your system's directory structure._

