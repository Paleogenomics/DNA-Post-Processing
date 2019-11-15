#!/bin/bash
# 31 Oct 2019
# A. Vershinina avershin@ucsc.edu https://github.com/avershinina
# Reference-guided aseembly of ancient mitochondrial genome
# Quickstart: ./mito_assembly_pipeline.sh path/to/sample.R1.fastq.gz path/to/sample.R2.fastq.gz samplename

###################
# README
###################

# This script complements the paper by Vershinina et al 201X "The case of an arctic wild ass highlights the utility of ancient DNA for validating problematic identifications in museum collections", currently under review at Molecular Ecology Resources.
# This script is a pipeline that utilizes a set of bioinformatic tools to process anicent DNA sequencing data with the aim of assembling a complete mitochondrial genome using a reference. This script works only for paired end Illumina data.

# To run the pipeline:
# Inside of the script, a user would need to include the path to reference genome, path to software, flags and settings. All these are hard-coded inside of the script. To run the pipeline, edit all corresponding variables manually using any script-friendly text editor, such as sublime-text or gedit. Make the script executable: chmod +x mito_assembly.sh

# To run the pipeline, input three variables: path to forward and reverse read files and arbitrary sample name (used to name output files). Run it as follows:

# ./mito_assembly_pipeline.sh path/to/sample.R1.fastq.gz path/to/sample.R2.fastq.gz samplename

# Steps of the pipeline:
# 1) Trim adapters and merge reads. 
# Note. Since we are processing ancient DNA, we assume that endogenious sequences are mainly short (30-70bp). This, however, may not be true for well preserved DNA. User should check fragment length distribution before deciding if most reads in a read pair are overlapping and whether they should be merged or not. Not all of the sequences are overlapping enough to merge them. Thus, the current version of the script utilizes both merged (short) and unmerged (long) reads. User may be interested in using only merged reads (cases of extremely degraded samples). If that is the case, the script should be manually modified according to user's needs.

# 2) Remove low complexity reads (such as stretches of AAAAs, ATATATA, etc).

# 3) Concatenate files: collect both merged and unmerged reads together.

# 4) Convert fastq into fasta, applying a phred quality filter.

# 5) Remove duplicated reads (with the aim to reduce % of PCR duplicates).

# 6) Run MIA: mapping-iterative-assembler using a reference mitochondrial genome.

###################
# Dependencies
###################

# SEQPREP https://github.com/jeizenga/SeqPrep2
# MIA https://github.com/mpieva/mapping-iterative-assembler/tree/5a7fb5afad735da7b8297381648049985c599874
# PRINSEQ http://prinseq.sourceforge.net/
# BBMAP https://github.com/BioInfoTools/BBMap
# FASTX_TOOLKIT http://hannonlab.cshl.edu/fastx_toolkit/index.html


###################
# Pipeline start
# (no need to edit)
###################

READ1=$1 	# Your data filenames and (full) file path. For example: /data/novaseq/AV000_L001_R1_001.fastq.gz. 
READ2=$2 	# As above: /data/novaseq/AV000_L001_R2_001.fastq.gz
SAMPLE=$3 	# arbitrary, used to name files

########################
# EDIT VARIABLES BELOW
########################

# Software (each is a path to a folder, not to an executable file)
SEQPREP=/path/to/SeqPrep2-master
MIA=/path/to/mia-1.0/src 	# Both MIA and MA can be installed from here https://github.com/mpieva/mapping-iterative-assembler/tree/5a7fb5afad735da7b8297381648049985c599874
PRINSEQ=/path/to/prinseq
BBMAP=/path/to/bbmap/sh
FASTX_TOOLKIT=/path/to/fastx_toolkit-0.0.13.2/src 	# http://hannonlab.cshl.edu/fastx_toolkit/index.html

# Settings - seqprep
SEQPREP_MIN_L=27 	# Removes unmappably short reads. 
SEQPREP_OVERLAP=15 	# Merge reads if they overlap by XXbp
SEQPREP_MIN_Q=15	# Low quality cutoff

# Settings - prinseq
COMPLEXITY_METHOD=dust			# can change to 'entropy'
COMPLEXITY_THRESHOLD=7			# see prinseq manual. Change if using 'entropy'.

# Settings - fatsx toolking
FX_Q=33 # relatively high cut-off

# Settings - MIA
ANCIENT_SUBMAT=/path/to/ancient.submat.txt # for details see here https://github.com/mpieva/mapping-iterative-assembler/blob/5a7fb5afad735da7b8297381648049985c599874/matrices/ancient.submat.txt
K=14 # k-mer size
REF=/path/to/mitogenomes/organism.fasta # reference mitochondria. For example: /mitogenomes/Equus_caballus_NC_001640 
REF_NAME=name_mito # arbitrary. For example eqcab_NC001640

###########################
# PIPELINE (no need to edit)
###########################

echo '---------------------------'
echo 'Starting the pipeline'
echo '---------------------------'

mkdir seqprep_output
mkdir mia_output

${SEQPREP}/SeqPrep2 -f ${READ1} -r ${READ2} -1 seqprep_output/${SAMPLE}_R1_unmerged.fastq.gz -2 seqprep_output/${SAMPLE}_R2_unmerged.fastq.gz -q ${SEQPREP_MIN_Q} -L ${SEQPREP_MIN_L} -A AGATCGGAAGAGCACACGTC -B AGATCGGAAGAGCGTCGTGT -s seqprep_output/${SAMPLE}_merged.fastq.gz -E seqprep_output/${SAMPLE}_readable_alignment.txt -o ${SEQPREP_OVERLAP} -d 1 -C ATCTCGTATGCCGTCTTCTGCTTG -D GATCTCGGTGGTCGCCGTATCATT >& seqprep_output/${SAMPLE}_SeqPrep_output.txt

gzip -dc seqprep_output/${SAMPLE}_merged.fastq.gz | perl ${PRINSEQ}/prinseq-lite.pl -fastq stdin -out_good stdout -out_bad null -lc_method ${COMPLEXITY_METHOD} -lc_threshold ${COMPLEXITY_THRESHOLD} -line_width 0 | gzip > seqprep_output/${SAMPLE}_merged.complexity_filtered.fastq.gz

${BBMAP}/reformat.sh in1=seqprep_output/${SAMPLE}_R1_unmerged.fastq.gz in2=seqprep_output/${SAMPLE}_R2_unmerged.fastq.gz out=seqprep_output/${SAMPLE}_unmerged_combined.fastq.gz

echo '---------------------------'
echo 'Seqprep and bbmap have finished running'
echo '---------------------------'

gzip -dc seqprep_output/${SAMPLE}_unmerged_combined.fastq.gz | perl ${PRINSEQ}/prinseq-lite.pl -fastq stdin -out_good seqprep_output/${SAMPLE}_unmerged_combined.complexity_filtered -out_bad null -lc_method ${COMPLEXITY_METHOD} -lc_threshold ${COMPLEXITY_THRESHOLD} -line_width 0

gzip seqprep_output/${SAMPLE}_unmerged_combined.complexity_filtered.fastq

# Concatenate merged and unmerged reads together
cat seqprep_output/${SAMPLE}_unmerged_combined.complexity_filtered.fastq.gz seqprep_output/${SAMPLE}_merged.complexity_filtered.fastq.gz > seqprep_output/${SAMPLE}_all_trimmed.complexity_filtered.fastq.gz

gzip -dc seqprep_output/${SAMPLE}_all_trimmed.complexity_filtered.fastq.gz | ${FASTX_TOOLKIT}/fastq_to_fasta/fastq_to_fasta -n -Q ${FX_Q} -r > seqprep_output/${SAMPLE}_all_trimmed.complexity_filtered.fasta
# collapse duplicates

perl ${PRINSEQ}/prinseq-lite.pl -fasta seqprep_output/${SAMPLE}_all_trimmed.complexity_filtered.fasta -out_good seqprep_output/${SAMPLE}_all_trimmed.complexity_filtered.dedup -out_bad null -derep 124 -line_width 0

echo '---------------------------'
echo "Prinseq and fastx have finished running"
echo 'Starting mito assembly'
echo '---------------------------'

${MIA}/mia -r ${REF} -f seqprep_output/${SAMPLE}_all_trimmed.complexity_filtered.dedup.fasta -c -C -U -s ${ANCIENT_SUBMAT} -i -F -k $K -m mia_output/${SAMPLE}.${REF_NAME}.maln

echo '---------------------------'
echo 'Mito assembly has finished running'
echo '---------------------------'

${MIA}/ma -M mia_output/${SAMPLE}.${REF_NAME}.maln.* -f 3 > mia_output/${SAMPLE}.${REF_NAME}.maln.F.mia_stats.txt
${MIA}/ma -M mia_output/${SAMPLE}.${REF_NAME}.maln.* -f 2 > mia_output/${SAMPLE}.${REF_NAME}.maln.F.mia_coverage_per_site.txt
${MIA}/ma -M mia_output/${SAMPLE}.${REF_NAME}.maln.* -f 5 > mia_output/${SAMPLE}.${REF_NAME}.maln.F.mia_consensus.fasta
${MIA}/ma -M mia_output/${SAMPLE}.${REF_NAME}.maln.* -f 41 > mia_output/${SAMPLE}.${REF_NAME}.maln.F.inputfornext.txt

gzip seqprep_output/*.fasta

echo "Everything has finished running"
echo "MIA assembly stats:"
echo ${SAMPLE}.${REF_NAME}.maln.F.mia_stats.txt
head -n 5 mia_output/${SAMPLE}.${REF_NAME}.maln.F.mia_stats.txt