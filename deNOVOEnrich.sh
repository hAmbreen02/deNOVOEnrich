#!/bin/bash

## deNOVOEnrich is a computational pipeline developed to efficiently detect and profile
## genome-wide somatic transposition insertions, as well as non-reference heritable 
## transposition events from targeted sequencing data generated through Transposon display sequencing (TEd-Seq).
## Contact: Heena Ambreen (h.ambreen@exeter.ac.uk)
## Genome organization Lab: https://genorglab.wordpress.com/ & Bousios Lab: https://www.sussex.ac.uk/lifesci/bousioslab/
## Department of Biosciences, University of Exeter & School of Life Sciences, University of Sussex


echo "
#################################################
#######         deNOVOEnrich              #######
#######                                   #######
#######   De novo Transposition events    #######
#######                                   #######
####### Somatic & Heritable TE insertions #######
#################################################
"
echo ""
echo "For queries, comments or suggestions, contact Heena Ambreen at h.ambreen@exeter.ac.uk"
echo ""
echo "============================================================================="
show_help(){
echo ""
echo "Usage: bash deNOVOEnrich.sh [Options]"
echo ""
echo "Required Options:"
echo " --genome			Full path to directory containing reference genome FASTA"
echo " --ref_TE			Full path to directory with TE-specific FASTA and BED files: "
echo "                               - {TEfam}_flanking_sequence.fa"
echo "                               - {TEfam}_flanking_sequence_2.fa"
echo "                               - {TEfam}_loci.bed"
echo " --rawRead1			Path to raw paired-end Read 1 FASTQ file (.fq)"
echo " --rawRead2			Path to raw paired-end Read 2 FASTQ file (.fq)"
echo " --Sample			Name of the sample "
echo " --adapter			Path to Illumina adapter sequences file (e.g., TruSeq3-PE.fa)"
echo " --outDir			Directory to store all output results"
echo " --TEfam			Transposon family of Interest" 
echo "				(Must match TE-specific Samplenames example: Use AtCopia93 for"
echo "				AtCopia93_loci.bed or AtCopia93_flanking_sequence.fa"
echo " --monomeric_repeat		Bed file with coordinates of monomeric repeats of reference genome"
echo " --CORES			Number of threads to use"
echo " --help |-h			Display this help message"
echo ""
echo "Example:"
echo "  bash deNOVOEnrich.sh --Sample col1 --TEfam AtCopia93 --genome /path/genome.fa --ref_TE /path/TE_files --rawRead1 A1_1.fq --rawRead2 A1_2.fq --adapter /path/TruSeq3-PE.fa --outDir ./results --CORES 8"

echo "================================================================================"
echo 
}

for arg in "$@"; do
    if [[ "$arg" == "--help" || "$arg" == "-h" ]]; then
        show_help
        exit 0
    fi
done


######    Software Pre-requisites      #######
######---------------------------------#######
# - cutadapt version 5.1                     #
# - fastp version 1.0.1	        	     #
# - seqtk version 1.3-r106                   #
# - trimmomatic version 0.39                 #
# - PEAR version 0.9.6			     #
# - bowtie2 version 2.4.4		     #
# - picard version 3.4.0                     #
# - samtools version 1.13                    #
# - bedtools version 2.30.0                  #
# - readtagger version 0.5.25                #
# Ensure all tools are installed and in PATH#
#####----------------------------------####### 


######  Arguments & Input parameters   #######
######---------------------------------#######
#  --genome  Full path to directory containing reference genome FASTA
#  --ref_TE [path]  Full path to directory with TE-specific FASTA and BED files
#  --rawRead1  Path to raw paired-end Read 1 FASTQ file (.fq)
#  --rawRead2  Path to raw paired-end Read 2 FASTQ file (.fq)
#  --Sample	[string]	Name of the sample
#  --adapter [file]	Path to Illumina adapter sequences file (e.g., TruSeq3-PE.fa)
#  --outDir         Output directory where results will be generated
#  --TEfam          Transposon family of Interest. The name should match the TE-specific filenames
#  --CORES          Number of threads to use 
#  --help | -h		Display this help message
#####----------------------------------####### 


######        Input Samples         #######
######----------------------------#######
# Genome : Reference genome in FASTA format 
#
# Raw paired-end sequencing Samples from TEd-Seq in FASTQ format
#
# {TEfam}_flanking_sequence.fa: 
# Reverse complement sequence from the primer site to near the 5'LTR extremity in FASTA format

# {TEfam}_flanking_sequence_2.fa: 
# Forward sequence from extremity of the 5'LTR to the primer site in FASTA format
#
# {TEfam}_loci.bed :
# Genomic coordinates of the native TE loci in the reference genome 
# Format: ($Chr $Start $End). Column $4 with names is optional

# Genome_monomeric_repeat.bed:
# Genomic coordinates of mono-nucleotide repeats in reference genome 
#####----------------------------------####### 


######		Default Values 		  #######
genome=""                     
ref_TE=""                     
rawRead1=""                   
rawRead2=""                   
Sample=""                     
adapter="/path/to/TruSeq3-PE.fa"   
outDir="./deNOVOEnrich_output"    
TEfam="AtCopiaX"
monomeric_repeat=""            
CORES=8 


#####     Check Argument        #####
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --Sample) Sample="$2"; shift ;;
        --CORES) CORES="$2"; shift ;;
        --genome) genome="$2"; shift ;;
        --rawRead1) rawRead1="$2"; shift ;;
        --rawRead2) rawRead2="$2"; shift ;;
        --ref_TE) ref_TE="$2"; shift ;;
        --adapter) adapter="$2"; shift ;;
        --outDir) outDir="$2"; shift ;;
        --TEfam) TEfam="$2"; shift ;;
        --monomeric_repeat) monomeric_repeat="$2"; shift ;;
        --help|-h) show_help; exit 0 ;;
        *) echo " Unknown parameter: $1"; show_help; exit 1 ;;
    esac
    shift
done


if [[ -z "$genome" || -z "$ref_TE" || -z "$rawRead1" || -z "$rawRead2" || -z "$Sample" || -z "$adapter" || -z "$outDir" || -z "$TEfam" || -z "$monomeric_repeat" || -z "$CORES" ]]; then
    echo " ERROR: One or more required arguments are missing."
    show_help
    exit 1
fi

if [ -d "$outDir" ]; then
    echo "  Output directory '$outDir' already exists. Overwriting..."
    rm -rf "$outDir"
fi
mkdir -p $outDir

#path to Trimmed reads directory
outTrim=$outDir/Trimmed
mkdir -p $outTrim


genome_name=$(basename "$genome")        
genome_base="${genome_name%.*}" 


echo "deNOVOEnrich is now running — analysing transposition events..."
echo ""
echo "============================================================================="

echo " Chosen options"
echo "Sample Name            : $Sample"
echo "Read 1 FASTQ           : $rawRead1"
echo "Read 2 FASTQ           : $rawRead2"
echo "Reference Genome       : $genome"
echo "TE Reference Directory : $ref_TE"
echo "TE Family              : $TEfam"
echo "Adapter File           : $adapter"
echo "Output Directory       : $outDir"
echo "Threads (CORES)        : $CORES"
echo "Monomeric Repeat BED   : $monomeric_repeat"
echo ""
echo "============================================================================="


##Run pipeline##
##Ia## Retrieving TE-specific reads and remove low Quality or read artefacts followed by Quality control and adapter trimming
echo ""
echo "Retrieving TE-specific reads and removing low Quality or artefactual reads from raw sequencing data"
echo ""
source activate denovoenrich

cutadapt --discard-untrimmed -g file:$ref_TE/${TEfam}_flanking_sequence.fa -e 0.2 -o $outTrim/sub_1.fastq -p $outTrim/sub_2.fastq $rawRead1 $rawRead2 -O 25 -j $CORES --nextseq-trim=20

cutadapt -A file:$ref_TE/${TEfam}_flanking_sequence_2.fa -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTGGCATCTA -e 0.1 -o $outTrim/sub_final_1.fastq -p $outTrim/sub_final_2.fastq $outTrim/sub_1.fastq $outTrim/sub_2.fastq -j $CORES --nextseq-trim=20 -m 20


fastp -i $outTrim/sub_final_1.fastq -I $outTrim/sub_final_2.fastq -o $outTrim/sub_final_1a.fastq -O $outTrim/sub_final_2a.fastq -g -x -w 7 -l 20 -y -A

read1_count=$(wc -l < "$outTrim/sub_final_1a.fastq")
read1_count=$((read1_count / 4))

read2_count=$(wc -l < "$outTrim/sub_final_2a.fastq")
read2_count=$((read2_count / 4))


##Ib## Fetch fastq-read ID from trimmed fastq {Sample}s obtained after fastp run

grep "@" $outTrim/sub_final_1a.fastq | sed 's/^@//; s/ /_/'g > $outTrim/sub_final_1a.fastqID.txt 
grep "@" $outTrim/sub_final_2a.fastq | sed 's/^@//; s/ /_/'g > $outTrim/sub_final_2a.fastqID.txt 

id1_count=$(wc -l < "$outTrim/sub_final_1a.fastqID.txt") ##number should match the one in sub_final_1a.fastq
id2_count=$(wc -l < "$outTrim/sub_final_2a.fastqID.txt") ##number should match the one in sub_final_2a.fastq

if [[ "$read1_count" -eq "$id1_count" && "$read2_count" -eq "$id2_count" ]]; then
    echo " Read count matches between trimmed FASTQ and ID lists."
else
    echo " Mismatch detected! Check FASTQ integrity or ID extraction."
    echo "Read1: $read1_count IDs vs $id1_count lines"
    echo "Read2: $read2_count IDs vs $id2_count lines"
fi



##Ic## Retrieve raw fastq reads

#zcat $rawRead1 | sed 's/[ \t]/_/g' > $outTrim/${Sample}_1_temp.fq
sed 's/[ \t]/_/g' $rawRead1 > $outTrim/${Sample}_1_temp.fq
gzip $outTrim/${Sample}_1_temp.fq


#zcat $rawRead2 | sed 's/[ \t]/_/g' > $outTrim/${Sample}_2_temp.fq
sed 's/[ \t]/_/g' $rawRead2 > $outTrim/${Sample}_2_temp.fq
gzip $outTrim/${Sample}_2_temp.fq


seqtk subseq $outTrim/${Sample}_1_temp.fq.gz $outTrim/sub_final_1a.fastqID.txt > $outDir/${Sample}_sub_final_1a_raw.fastq
seqtk subseq $outTrim/${Sample}_2_temp.fq.gz $outTrim/sub_final_2a.fastqID.txt > $outDir/${Sample}_sub_final_2a_raw.fastq




##Id## Quality control and adapter trimming

trimmomatic PE -threads $CORES -trimlog $outDir/${Sample}_raw_sub_final_${TEfam}_norm1.trimlog -summary $outDir/${Sample}_raw_sub_final_${TEfam}_norm.summary $outDir/${Sample}_sub_final_1a_raw.fastq $outDir/${Sample}_sub_final_2a_raw.fastq -baseout $outDir/${Sample}_raw_sub_final_${TEfam}_norm ILLUMINACLIP:${adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

echo ""
echo "TE-specific high quality reads retrieved"
echo ""
echo "============================================================================="


##IIa## Merge PE reads into single reads

echo ""
echo "Generating long contigs"
echo ""
pear -f $outDir/${Sample}_raw_sub_final_${TEfam}_norm_1P -r $outDir/${Sample}_raw_sub_final_${TEfam}_norm_2P -o $outDir/${Sample}_raw_sub_final_${TEfam}_norm -v 5 -m 300 -n 50 -j $CORES -y 100M


echo "Paired-end reads merged"
echo ""
echo "============================================================================="
echo ""

##IIIa## Mapping on reference

echo ""
echo "Mapping on reference genome"
echo ""
          
bowtie2-build --threads $CORES $genome $outDir/$genome_base

bowtie2 -x $outDir/$genome_base -U $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.fastq -S $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sam --local --very-sensitive --threads $CORES

samtools view -bhS -@ $CORES $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sam | samtools sort -@ $CORES -o $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.bam

echo "Reference mapping done"
echo ""

##IIIb## Retrieving unique and clipped alignments (minimum >40nt clipping)

samtools view -H --threads $CORES $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.bam > $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.uniq.bam

samtools view -@ $CORES -F 4 $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.bam | awk -v MR=0.8 'BEGIN{FS=OFS="\t"} {for(i=12;i<=100;i++){if($i~/AS:i:/){as=substr($i,6)};if($i~/XS:i:/){xs=substr($i,6);break}};if(xs/1<as*MR){print $0}}' >> $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.uniq.bam

samtools view -H --threads $CORES $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.uniq.bam > $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.uniq.tmp.header

awk '$6~/^[4-8][0-9]S/ || $6~/[4-8][0-9]S$/ || $6~/^1[0-9][0-9]S/ || $6~/1[0-9][0-9]S$/' $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.uniq.bam | cat $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.uniq.tmp.header - | samtools view -@ 40 -bh - | samtools sort -@ 40 -o $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.uniq.split.bam -

echo "You now have Unique and Clipped ALignemnts from reference mapping"
echo ""
echo "============================================================================="
echo ""
##IIIc##Marking duplicates on unique and clipped alignment Sample


picard MarkDuplicates REMOVE_DUPLICATES=false I=$outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.uniq.split.bam O=$outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.uniq.split.dedup.bam M=$outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.uniq.split.dedup.metrics

samtools index $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.uniq.split.dedup.bam $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.uniq.split.dedup.bai

echo ""
echo "Part1: Reference Samples generated, Moving on to TE mapping......."
echo ""
echo "============================================================================="

##IVa## Mapping the QC filtered reads on TE extremity and marking duplicates

echo ""
echo "Mapping on reference TE"

bowtie2-build --threads $CORES $ref_TE/${TEfam}_flanking_sequence_2.fa $outDir/${TEfam}_extremity

bowtie2 -x $outDir//${TEfam}_extremity -U $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.fastq -S $outDir/${Sample}_raw_sub_final_norm.assembled_${TEfam}_ext.sam --local --very-sensitive --threads $CORES

samtools view -bhS -@ $CORES $outDir/${Sample}_raw_sub_final_norm.assembled_${TEfam}_ext.sam | samtools sort -@ $CORES -o $outDir/${Sample}_raw_sub_final_norm.assembled_${TEfam}_ext.sorted.bam


picard MarkDuplicates REMOVE_DUPLICATES=false I=$outDir/${Sample}_raw_sub_final_norm.assembled_${TEfam}_ext.sorted.bam O=$outDir/${Sample}_raw_sub_final_norm.assembled_${TEfam}_ext.sorted.dedup.bam M=$outDir/${Sample}_raw_sub_final_norm.assembled_${TEfam}_ext.sorted.dedup.metrics

samtools index $outDir/${Sample}_raw_sub_final_norm.assembled_${TEfam}_ext.sorted.dedup.bam $outDir/${Sample}_raw_sub_final_norm.assembled_${TEfam}_ext.sorted.dedup.bai

echo ""
echo "You now have Reference and TE ALignments"
echo ""
echo "============================================================================="

#conda deactivate
##Va## Tagging Alignments with Readtagger (target: Reference; source: TE_extremity)

echo ""
echo "Tagging Genome:TE alignments"

source activate readtagger
readtagger -t $outDir/${Sample}_raw_sub_final_${TEfam}_norm.assembled.sorted.uniq.split.dedup.bam -s $outDir/${Sample}_raw_sub_final_norm.assembled_${TEfam}_ext.sorted.dedup.bam -o $outDir/${Sample}_${TEfam}.uniq.split.readtagged.bam --cores $CORES 

#conda deactivate

##Vb## Retrieving read-tagged alignment records for peak calling

source activate denovoenrich
samtools view -H --threads $CORES $outDir/${Sample}_${TEfam}.uniq.split.readtagged.bam > $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bam

samtools view $outDir/${Sample}_${TEfam}.uniq.split.readtagged.bam | grep "${TEfam}" >> $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bam



echo ""
echo "============================================================================="
echo "Calling non-reference somatic and heritable TE insertions"
echo ""


##VIa## Peak calling
##Using the bamtobed Samples with PCR duplicates: this will allow us to retrive all the unique reads which will be putative insertions and if a particular insertion events has clonal duplications or not. One with very high pcr duplicates will be fixed insertion. Other with low or no duplicates will be somatic insertions

bedtools bamtobed -i $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bam -split > $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.bed

sort -k1,1 -k2,2n $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.bed  | awk ' $1 !="ChrM" && $1 != "ChrC" {print $0}' > $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.sorted.bed ## 6 column bed: chr start stop readID score strand

awk '{print $1"\t"$2"\t"$3"\t"$6}' $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.sorted.bed | awk '{count[$0]++} END {for (line in count) print line, count[line]}'| sort -k1,1 -k2,2n | sed 's/ /\t/g' > $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.sorted.withPCRdupCount.bed 

awk '$4=="+" {print $1"\t"$2"\t"$4"\t"$5}' $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.sorted.withPCRdupCount.bed | sort -k1,1 -k2,2n | bedtools groupby -i stdin -g 1,2,3 -c 4 -o sum,count,collapse > $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.startsite.bed

awk '$4=="-" {print $1"\t"$3"\t"$4"\t"$5}' $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.sorted.withPCRdupCount.bed | sort -k1,1 -k2,2n | bedtools groupby -i stdin -g 1,2,3 -c 4 -o sum,count,collapse >> $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.startsite.bed

sort -k1,1 -k2,2n $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.startsite.bed > $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.startsite.sorted.bed ##this will generate a 6 bed Sample for start positions: chr start_position strand sumofallreads_includePCRduplicates countofreadswithcommonstartbutdiffstopsite comma_sperated_list_of_ind_pcrdup_counts




##VIb## Calling fixed insertions: conditional filtering for non-reference fixed loci: threshold coverage >=30, minimum number of unique reads at the insertion site >=5 and from the unique reads atleast 3 of them will have pcr duplicates of 3 or more #


 
awk -F'\t' '{split($6, values, ","); count = 0; for (i in values) {if (values[i] >= 3) {count++;}} if ($4 >= 30 && $5 >= 5 && count >= 3) {print $0;}}' $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.startsite.sorted.bed > $outDir/${Sample}_${TEfam}.putative.fixed_insertions.bed

##convert into bed first here you have only start coordinates also you need to take care of the strand info in the last call whether the TE seq is direct or reverse complement

awk '{{if ($3=="-") print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6} if( $3=="+")  print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6}' $outDir/${Sample}_${TEfam}.putative.fixed_insertions.bed > $outDir/${Sample}_${TEfam}.putative.fixed_insertions.temp1.bed #7 column bed: chr start stop  strand total_coverage Total_unique reads duplicates_distribution_per_unique_read 


##VIc## Filtering putative insertions: removing insertions within 1kb of native loci and those within 100bp of mono-nucleotide trails


bedtools window -v -w 1000 -a $outDir/${Sample}_${TEfam}.putative.fixed_insertions.temp1.bed  -b $ref_TE/${TEfam}_loci.bed > $outDir/${Sample}_${TEfam}.fixed_insertions.temp2.bed

bedtools window -v -w 50 -a $outDir/${Sample}_${TEfam}.fixed_insertions.temp2.bed -b $monomeric_repeat > $outDir/${Sample}_${TEfam}.fixed_insertions.final.bed ##use this file in case bed format is required for intersection with other features

awk '{{if ($4=="-") print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} if( $4=="+")  print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' $outDir/${Sample}_${TEfam}.fixed_insertions.final.bed > $outDir/${Sample}_${TEfam}.fixed_insertions.final.startsite.bed ##use this Sample if a specific insertion site is required

heritable_count=$(wc -l < "$outDir/${Sample}_${TEfam}.fixed_insertions.final.bed")

echo ""
echo ""$heritable_count putative heritable/fixed insertions of $TEfam detected for $Sample""
echo ""
echo "============================================================================="
echo ""

##VIIa##Calling somatic insertions: conditional filtering for somatic insertions

## rare somatic insertions where no pcr duplicates are found and a single read represent the insertion site
awk '$5==1 && $4==1 {print $0}' $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.startsite.sorted.bed > $outDir/${Sample}_${TEfam}.putative_rare_soma_insertions.bed
awk '{{if ($3=="-") print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6} if( $3=="+")  print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6}' $outDir/${Sample}_${TEfam}.putative_rare_soma_insertions.bed > $outDir/${Sample}_${TEfam}.putative_rare_soma_insertions.temp1.bed
bedtools window -v -w 1000 -a $outDir/${Sample}_${TEfam}.putative_rare_soma_insertions.temp1.bed  -b $ref_TE/${TEfam}_loci.bed > $outDir/${Sample}_${TEfam}.putative_rare_soma_insertions.temp2.bed
bedtools window -v -w 50 -a $outDir/${Sample}_${TEfam}.putative_rare_soma_insertions.temp2.bed -b $monomeric_repeat > $outDir/${Sample}_${TEfam}_rare_soma_insertions.final.bed ##use this Sample in case bed format is required for intersection with other features
awk '{{if ($4=="-") print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} if( $4=="+")  print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' $outDir/${Sample}_${TEfam}_rare_soma_insertions.final.bed > $outDir/${Sample}_${TEfam}_rare_soma_insertions.final.startsite.bed ##use this Sample if a specific insertion site is required


##recent somatic insertions where at least 1 or more PCR duplicate is available for the insertion site
awk '$5==1 && $4>1 {print $0}' $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.startsite.sorted.bed > $outDir/${Sample}_${TEfam}.putative_recent_soma_insertions.bed
awk '{{if ($3=="-") print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6} if( $3=="+")  print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6}' $outDir/${Sample}_${TEfam}.putative_recent_soma_insertions.bed > $outDir/${Sample}_${TEfam}.putative_recent_soma_insertions.temp1.bed
bedtools window -v -w 1000 -a $outDir/${Sample}_${TEfam}.putative_recent_soma_insertions.temp1.bed  -b $ref_TE/${TEfam}_loci.bed > $outDir/${Sample}_${TEfam}.putative_recent_soma_insertions.temp2.bed
bedtools window -v -w 50 -a $outDir/${Sample}_${TEfam}.putative_recent_soma_insertions.temp2.bed -b $monomeric_repeat > $outDir/${Sample}_${TEfam}_recent_soma_insertions.final.bed ##use this Sample in case bed format is required for intersection with other features
awk '{{if ($4=="-") print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} if( $4=="+")  print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' $outDir/${Sample}_${TEfam}_recent_soma_insertions.final.bed > $outDir/${Sample}_${TEfam}_recent_soma_insertions.final.startsite.bed ##use this Sample if a specific insertion site is required


##primary somatic insertions (occurred very early in development/experimental phase) where more than 1 unique reads with 1 or more pcr duplicates are available 
awk '$5 < 5 && $5 > 1 && $4 >= 2 * $5 && $4 < 30 {print $0}' $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bamtobed.startsite.sorted.bed > $outDir/${Sample}_${TEfam}.putative_early_soma_insertions.bed
awk '{{if ($3=="-") print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6} if( $3=="+")  print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6}' $outDir/${Sample}_${TEfam}.putative_early_soma_insertions.bed > $outDir/${Sample}_${TEfam}.putative_early_soma_insertions.temp1.bed
bedtools window -v -w 1000 -a $outDir/${Sample}_${TEfam}.putative_early_soma_insertions.temp1.bed  -b $ref_TE/${TEfam}_loci.bed > $outDir/${Sample}_${TEfam}.putative_early_soma_insertions.temp2.bed
bedtools window -v -w 50 -a $outDir/${Sample}_${TEfam}.putative_early_soma_insertions.temp2.bed -b $monomeric_repeat > $outDir/${Sample}_${TEfam}_early_soma_insertions.final.bed ##use this Sample in case bed format is required for intersection with other features
awk '{{if ($4=="-") print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} if( $4=="+")  print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' $outDir/${Sample}_${TEfam}_early_soma_insertions.final.bed > $outDir/${Sample}_${TEfam}_early_soma_insertions.final.startsite.bed ##use this Sample if a specific insertion site is required

##VIIb##combine all soma insertions: Total somatic insertions

cat $outDir/${Sample}_${TEfam}_*_soma_insertions.final.bed | sort -k1,1 -k2,2n > $outDir/${Sample}_${TEfam}_All_soma_insertions.final.bed

awk '$5<=5 {print $0}' $outDir/${Sample}_${TEfam}_All_soma_insertions.final.bed > $outDir/${Sample}_${TEfam}.somatic_insertions.final.bed

somatic_count=$(wc -l < "$outDir/${Sample}_${TEfam}.somatic_insertions.final.bed")

echo ""
echo ""$somatic_count putative somatic insertions of $TEfam detected for $Sample""
echo ""
echo "============================================================================="
echo ""

samtools view -hb $outDir/${Sample}_${TEfam}.uniq.split.readtagged.final.bam > $outDir/${Sample}_${TEfam}.uniq.split.tagged_alignments.final.bam

samtools index $outDir/${Sample}_${TEfam}.uniq.split.tagged_alignments.final.bam


#conda deactivate 
rm $outDir/*.temp*.bed
rm -r $outTrim

temp=$outDir/tempDir
mkdir -p $temp
mv $outDir/* $temp

finalDir=$outDir/Final_outputs
mkdir -p $finalDir
mv $temp/${Sample}_${TEfam}.somatic_insertions.final.bed $temp/${Sample}_${TEfam}.fixed_insertions.final.bed $temp/${Sample}_${TEfam}.uniq.split.tagged_alignments.final.bam $temp/${Sample}_${TEfam}.uniq.split.tagged_alignments.final.bam.bai $outDir/Final_outputs/ 

echo "============================================================================="
echo ""
echo "Analysis for TE insertion calling with deNOVOEnrich finished"
echo ""
echo ""
echo "Thank you for using the tool!"
echo ""
echo "Happy TEing!!!"
echo ""
echo "============================================================================="
echo ""
