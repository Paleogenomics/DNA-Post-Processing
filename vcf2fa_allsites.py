#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
import argparse
import gzip
import subprocess
import random
from collections import defaultdict
import os
"""
vcf2fa_allsites.py
Created by Nathan Schaefer on 05/25/16 at 11:58

Given a VCF file and a FASTA index file for the reference genome, 
outputs 2 'haplotypes' for each individual in the VCF file.

If there are heterozygous sites, alleles are randomly assigned to one or
    the other haplotype. In other words, this does not phase data, and neither
    'haplotype' output is a true haplotype.

VCF files must be created using GATK's "EMIT_ALL_SITES" options. In other words,
    if a site is missing from the file, it is assumed that that site failed
    to pass quality filters -- NOT that it is a reference allele. Any sites missing
    from the VCF will appear as ambiguous/"N" bases in the output files.
    
Requires that SAMTools and BGZIP are installed and in your PATH
    (to check, try which samtools and which bgzip)

Also expects uncompressed VCF input via stdin. If you want to do filtering
    before this, pipe from vcftools (i.e. vcftools --gzvcf $FILE --remove-filtered-all
    --minQ 50 --recode --stdout | $thisprogram)
    Or grab a desired region via tabix (i.e. tabix $FILE scaffold1:0-1000000 | $thisprogram)
    Or to grab the whole file (i.e. zcat $FILE | $thisprogram)
"""

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    
    out_opts = parser.add_mutually_exclusive_group(required=True)
    
    out_opts.add_argument("--recode", "-r", help=\
        "Output VCF format",\
        action="store_true",\
        required=False)
        
    fa_opts = out_opts.add_argument_group()
    fa_opts.add_argument("--fai", "-f", help=\
        "The path to the FAI file for the reference genome",
        type=parse_fai)
    fa_opts.add_argument("--output_prefix", "-o", help=\
        "The prefix of files to create (will be followed by sample name_hapX.fa)")
    
    parser.add_argument("--vq", "-v", help=\
        "Minimum variant quality to allow through filter.",
        type=float,
        default=0)
    parser.add_argument("--gq", "-g", help=\
        "Minimum genotype quality to allow through filter.",
        type=float,
        default=0)
    parser.add_argument("--mq", "-m", help=\
        "Minimum root-mean-squared map quality to allow through filter.",
        type=float,
        default=0)
    parser.add_argument("--mindepth", help=\
        "Minimum depth that ALL individuals in the file must have to pass filter.",
        type=int,
        default=0)
    parser.add_argument("--maxdepth", help=\
        "Maximum depth below which ALL individuals in the file must fall to pass filter.",
        type=int,
        default=99999)
    parser.add_argument("--depth", "-d", help=\
        "Provide individual depth filters for every individual in the file as a string. "
        "String should be in the format INDV=MIN-MAX (comma separated) and should have an entry for every "
        "individual in the file (or it will fail).",
        required=False,
        type=parse_depthstr)
    parser.add_argument("--no_missing", "-nm", help=\
        "This flag will cause ALL individuals to have N in their sequences at sites "
        "where ANY individual has missing (or filtered) genotype data. The default "
        "is to only put N into sequences for those individuals with missing data.",
        action="store_true",
        default=False)
    
    parsed = parser.parse_args()
    if not parsed.recode:
        # Make sure both pieces of information needed are present for writing
        # FASTA files.
        if options.output_prefix is None:
            print("ERROR: you must provide an output prefix for writing FASTA files.", file=sys.stderr)
            exit(1)
        if options.fai is None:
            print("ERROR: you must provide the path to an index (.fai) of the reference FASTA file in order to output FASTA files.", file=sys.stderr)
            exit(1)
    return parsed

def parse_depthstr(depths):
    items = depths.split(',')
    # Store depth as dict. Keys are individual names, values are tuples of
    # (mindepth, maxdepth)
    depth = {}
    for item in items:
        try:
            indv, minmax = item.split('=')9,
            mincov, maxcov = minmax.split('-')
            depth[indv] = (int(mincov), int(maxcov))
        except:
            print("ERROR: the depth string you have provided is invalid.",file=sys.stderr)
            print("It should be in the form INDV=MIN-MAX for each individual (comma separated).",file=sys.stderr)
            print("MIN and MAX must be integer values.", file=sys.stderr)
            exit(1)
    return depth
    
def parse_fai(filename):
    chrlens = {}
    f = open(filename, 'r')
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        chrlens[data[0]] = int(data[1])
    f.close()
    return chrlens
    
def main(args):
    """Main method"""
    options = parse_args()

    samples = {}
    
    # This will store the sequence of the current chromosome/scaffold, as a 
    # list (so it's mutable), keyed to individual's sample names (not indices)
    cur_seqs_1 = defaultdict(list)
    cur_seqs_2 = defaultdict(list)
    
    # Store files for writing out haplotypes
    files_1 = {}
    files_2 = {}
    
    # Store file names so we can gzip them.
    filenames = []
    
    prev_chrom = None
    
    sites_tot = 0
    sites_valid = 0
    
    for line in sys.stdin:
        line = line.rstrip()
        if len(line) < 2:
            continue
        elif line[0:2] == "##":
            if options.recode:
                print(line)
            continue
        elif line[0] == "#":
            if options.recode:
                print(line)
            # Header line.
            data = line.split('\t')
            for sampleIndex in range(9, len(data)):
                sample_name = data[sampleIndex]
                if options.depth is not None and sample_name not in options.depth:
                    print("ERROR: you did not provide depth cutoffs for individual {}".format(sample_name), file=sys.stderr)
                    exit(1)
                samples[sampleIndex] = sample_name

                # Determine file names
                if not options.recode:
                    if os.path.isdir(options.output_prefix):
                        if options.output_prefix[-1] == "/":
                            options.output_prefix = options.output_prefix[0:-1]
                        fn1 = "{}/{}_hap1.fa".format(options.output_prefix, sample_name)
                        fn2 = "{}/{}_hap2.fa".format(options.output_prefix, sample_name)
                    else:
                        fn1 = "{}_{}_hap1.fa".format(options.output_prefix, sample_name)
                        fn2 = "{}_{}_hap2.fa".format(options.output_prefix, sample_name)
                    
                    filenames.append(fn1)
                    filenames.append(fn2)
                    
                    # Create files                    
                    files_1[sample_name] = open(fn1, 'w')
                    files_2[sample_name] = open(fn2, 'w')
                 
        else:
            sites_tot += 1
            
            data = line.split('\t')
            chrom = data[0]
            if chrom != prev_chrom:
                if prev_chrom is not None and not options.recode:
                    # Print out data.
                    for sample_name in cur_seqs_1:
                        print(">{}".format(prev_chrom), file=files_1[sample_name])
                        print("".join(cur_seqs_1[sample_name]), file=files_1[sample_name])
                        print(">{}".format(prev_chrom), file=files_2[sample_name])
                        print("".join(cur_seqs_2[sample_name]), file=files_2[sample_name])
                
                if not options.recode:
                    # Create new empty sequences.
                    # Each sequence is a list of Ns the same length as the 
                    # chromosome/scaffold, according to the FASTA index.
                    for indv_index in samples:
                        cur_seqs_1[samples[indv_index]] = ["N"] * options.fai[chrom]
                        cur_seqs_2[samples[indv_index]] = ["N"] * options.fai[chrom]
                    
                prev_chrom = chrom
            
            # Start a sequence if we haven't already.
            pos = int(data[1])
            # Convert to 0-based coordinates
            pos -= 1
            
            ref = data[3]
            alt = data[4].split(',')
            
            if data[5] == '.':
                qual = -1
            else:
                qual = float(data[5])
            
            filter_string = data[6]
            
            # At this point, we can skip the site if there is no alt allele,
            # variant quality is low, the built-in filter did not pass, or
            # any alt allele is an indel.
            if filter_string != "." and filter_string != "PASS":
                continue
            if qual < options.vq:
                continue
            for alt_allele in alt:
                if alt_allele == "." or len(alt_allele) > 1:
                    continue
            
            mq = None
            info = data[7].split(';')
            for item in info:
                if item != '.':
                    key, val = item.split('=')
                    if key == "MQ":
                        mq = float(val)
                        break
            
            # Check to see if we passed root-mean-squared map quality filter.
            if mq is not None and options.mq is not None and mq < options.mq:
                continue
            
            gt_field = None
            dp_field = None
            gq_field = None
            
            fields = data[8].split(':')
            for index, fld in enumerate(fields):
                if fld == "GT":
                    gt_field = index
                elif fld == "DP":
                    dp_field = index
                elif fld == "GQ":
                    gq_field = index
            
            # At this point, we can skip the site if we are missing genotype,
            # depth, or genotype quality information.
            if gt_field is None or dp_field is None:
                continue
            
            # Now we can get individual genotypes.
            missing_genotype = False
            
            alleles = [ref] + alt
            
            sample_alleles = {}
            
            for indv_index in samples:
                
                indv_name = samples[indv_index]
                
                if options.recode:
                    # Fill in missing data first; update later if we find non-missing
                    # data.
                    sample_alleles[indv_index] = ('.', '.')
                    
                sampledata = data[indv_index].split(':')
                if len(sampledata) < max(gt_field, dp_field, gq_field) + 1:
                    # Missing data.
                    missing_genotype = True
                    continue

                gt = sampledata[gt_field]
                dp = sampledata[dp_field]
                if dp == '.':
                    dp = 0
                else:
                    dp = int(dp)
                    
                if gq_field is not None:
                    gq = sampledata[gq_field]
                    if gq == '.':
                        gq = 0
                    else:
                        gq = int(gq)
                elif alt[0] == ".":
                    # If it's a homozygous reference site, GATK doesn't seem
                    # capable of computing genotype qualities, so we won't 
                    # require them (fake the program into thinking we passed
                    # filter)
                    gq = options.gq
                else:
                    gq = -1
                    
                if '/' in gt:
                    gt = gt.split('/')
                else:
                    gt = gt.split('|')
                
                
                # Check if the individual meets all requirements; otherwise, 
                # leave "N"
                if gq >= options.gq and dp >= options.mindepth and \
                    dp <= options.maxdepth and gt[0] != '.' and gt[1] != '.':
                    
                    gt[0] = int(gt[0])
                    gt[1] = int(gt[1])
                
                    # Also check individual coverage depth filter, if given
                    if options.depth is not None and (dp < options.depth[indv_name][0] \
                        or dp > options.depth[indv_name][1]):
                        # Individual coverage filter failed.
                        missing_genotype = True
                    else:
                        # Store alleles.
                        if options.recode:
                            sample_alleles[indv_index] = (gt[0], gt[1])
                        else:
                            # Randomly assign heterozygous sites to one or the other
                            # haplotype
                            if gt[0] != gt[1]:
                                if random.random() < 0.5:
                                    sample_alleles[indv_name] = (alleles[gt[0]], alleles[gt[1]])
                                else:
                                    sample_alleles[indv_name] = (alleles[gt[1]], alleles[gt[0]])
                            else:
                                sample_alleles[indv_name] = (alleles[gt[0]], alleles[gt[1]])
                else:
                    missing_genotype = True
            
            if options.recode:
                
                # At this point, we've passed all site filters. We can now 
                # stick the genotypes back into "data" and print.
                for indv_index in samples:
                    # Rebuild based on "INFO" tag.
                    indv_info = data[indv_index].split(':')
                    indv_info[gt_field] = '{}{}{}'.format(sample_alleles[indv_index][0],\
                        indv_info[gt_field][1], sample_alleles[indv_index][1])
                    data[indv_index] = ':'.join(indv_info)
            
                print("\t".join(data))
            
            # If we don't require every individual to have good/usable data,
            # go ahead and add all stored variants to the sequences.
            if not (options.no_missing and missing_genotype):
                if options.recode:
                    pass
                else:
                    sites_valid += 1
                    for indv_name in sample_alleles:
                        cur_seqs_1[indv_name][pos] = sample_alleles[indv_name][0]
                        cur_seqs_2[indv_name][pos] = sample_alleles[indv_name][1]
            
            #print(sites_valid/sites_tot)
    
    if not options.recode:
        # Print out data for the last chromosome/scaffold.
        if prev_chrom is not None:
            for sample_name in cur_seqs_1:
                print(">{}".format(prev_chrom), file=files_1[sample_name])
                print("".join(cur_seqs_1[sample_name]), file=files_1[sample_name])
                print(">{}".format(prev_chrom), file=files_2[sample_name])
                print("".join(cur_seqs_2[sample_name]), file=files_2[sample_name])
        
        # Close all files.
        for sample_name in files_1:
            files_1[sample_name].close()
            files_2[sample_name].close()
        
        # BGZip and samtools faidx all FASTA files.
        for filename in filenames:
            subprocess.call(['bgzip', filename])
            subprocess.call(['samtools', 'faidx', '{}.gz'.format(filename)])
                
if __name__ == "__main__" :
    sys.exit(main(sys.argv))

