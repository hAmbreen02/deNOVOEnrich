# deNOVOEnrich
Identification of _de novo_ transposon insertion sites from targeted sequencing data generated through Transposon display Sequencing (TEd-Seq)

<img src="images/denonvoEnrich_logo.jpg" alt="denovoEnrich logo" width="300" height="300"/>

## Overview

deNOVOEnrich is a computational pipeline developed to efficiently detect and profile genome-wide somatic transposition insertions, as well as non-reference heritable transposition events from targeted sequencing data generated through Transposon display sequencing (TEd-Seq). 

The pipeline leverages high-confidence split read alignments to identify TE:genome break sites, enabling accurate detection of bona-fide non-reference somatic and heritable new insertions of a transposon. Stringent filters were applied at several critical steps of the pipeline to ensure removal of false positives and reduction of background noise signal.

The details of the pipeline have been described in Ambreen et al. 

## Requirements

- Linux 86X64 Systems
- cutadapt (https://cutadapt.readthedocs.io/en/stable/)
- fastp (https://github.com/OpenGene/fastp)
- seqtk (https://github.com/lh3/seqtk)
- trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)
- PEAR (https://github.com/tseemann/PEAR)
- bowtie2 (https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- picard (https://broadinstitute.github.io/picard/)
- samtools (https://www.htslib.org/)
- bedtools (https://bedtools.readthedocs.io/en/latest/index.html)
- readtagger (https://pypi.org/project/readtagger/)

Ensure that all required dependencies are properly installed in the working environment to avoid errors during run.






