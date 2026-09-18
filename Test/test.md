# Running deNOVOEnrich on Test samples 
This directory contains test files to validate and demonstrate the usage of the deNOVOEnrich pipeline.

## Note on raw sequencing input files
The paired-end test FASTQ files (Test_1.fq.bz2 & Test_2.fq.bz2) are compressed with bzip2. Before running the pipeline, these files should be extracted, as the pipeline **does not accept .bz2 files** directly.

## Reference genome
The test run uses chromosome 1 of the Arabidopsis thaliana **Col-CEN reference** assembly.

## Mononucleotide file
The BED file, **Arabidopsis_mononucl_repeat.bed**, includes genomic coordinates for mononucleotide repeats for chromosome 1 of the Col-CEN reference.

# Example command to run the test files

<pre> bash deNOVOEnrich.sh --Sample Test --TEfam AtCopia93 --genome /path/genome.fa --ref_TE /path/TE_files --rawRead1 Test_1.fq --rawRead2 Test_2.fq --adapter /path/TruSeq3-PE.fa --somatic 5 --heritable 30 --outDir ./results --mononucl_repeat /path/mononucl_repeats.bed --CORES 8 </pre>

_Please ensure you update all file paths to match your system's directory structure._

