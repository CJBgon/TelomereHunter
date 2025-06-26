#!/usr/bin/python

# Copyright 2015 Lina Sieverling, Philip Ginsbach, Lars Feuerbach

# This file is part of TelomereHunter.

# TelomereHunter is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# TelomereHunter is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with TelomereHunter.  If not, see <http://www.gnu.org/licenses/>.


import os
import sys
import subprocess
import re
import pysam

#################################################################################
### filters telomere reads from an input BAM or CRAM file                    ###
### counts the total number of reads mapped to each chromosome band          ###
### determines the gc content distribution of the reads in the input file    ###
#################################################################################

def filter_telomere_reads(bam_file, band_file, out_dir, pid, sample, repeat_threshold_calc, repeat_threshold_set, mapq_threshold, repeats, consecutive_flag, remove_duplicates, input_reference_fasta=None):

    ################################################
    ### get patterns and make regular expression ###
    ################################################

    patterns_regex_forward = ""
    patterns_regex_reverse = ""

    for repeat in repeats:
        patterns_regex_forward += repeat + "|"
        patterns_regex_reverse += getReverseComplement(repeat) + "|"

    patterns_regex_forward = patterns_regex_forward[:-1]
    patterns_regex_reverse = patterns_regex_reverse[:-1]

    #########################
    ### open file handles ###
    #########################

    if bam_file.endswith(".bam"):
        file_mode = "rb"
        bamfile = pysam.Samfile(bam_file, file_mode)
    elif bam_file.endswith(".cram"):
        file_mode = "rc"
        if input_reference_fasta is None:
            raise ValueError("A reference FASTA file must be provided for CRAM input.")
        bamfile = pysam.Samfile(bam_file, file_mode, reference_filename=input_reference_fasta)
    else:
        raise ValueError("Input file must be .bam or .cram")

    filtered_file_path = out_dir + "/" + pid + "_filtered.bam"
    filtered_file = pysam.Samfile(filtered_file_path, "wb", template=bamfile)

    ######################################################
    ### check if repeat threshold has already been set ###
    ######################################################

    if repeat_threshold_calc != 'n':
        repeat_threshold = repeat_threshold_calc

    ############################
    ### make chromosome list ###
    ############################

    references = bamfile.references

    if references[0][0:3] == "chr":
        bam_chr_prefix = "chr"
    else:
        bam_chr_prefix = ''

    chromosome_list = [str(i) for i in range(1, 23)] + ["X", "Y"]
    chromosome_list_with_prefix = [bam_chr_prefix + chr for chr in chromosome_list]

    ################################################
    ### make band, spectrum and gc content lists ###
    ################################################

    bands_list = {chr: {"band_name": [], "end": []} for chr in chromosome_list + ["unmapped"]}
    spectrum_list = {chr: {} for chr in chromosome_list}

    for line in open(band_file, "r"):
        try:
            line = line.rstrip().split()
            end = line[2]
            band_name = line[3]

            if line[0][:3] == "chr":
                chrom_name = line[0][3:]
            else:
                chrom_name = line[0]

            bands_list[chrom_name]["band_name"].append(band_name)
            bands_list[chrom_name]["end"].append(int(end))
            spectrum_list[chrom_name][band_name] = {"reads": 0}

        except:
            print("Invalid line in banding file: \'" + " ".join(line) + "\'")

    spectrum_list["unmapped"] = {"unmapped": {"reads": 0}}
    bands_list["unmapped"]["band_name"].append("unmapped")
    bands_list["unmapped"]["end"].append(0)

    gc_content_list = {gc_content: 0 for gc_content in range(0, 101)}

    #############################
    ### loop through file ###
    #############################

    lastChromosome = ''
    chromosomeLsEnd = None
    chromosomeLsBand = None
    spectrumTemp = None
    spectrumTemp2 = None
    spectrumUnmapped = spectrum_list["unmapped"]["unmapped"]
    sequence = ""
    chr_offset = len(bam_chr_prefix)
    i = 0
    band = ''

    for read in bamfile.fetch(until_eof=True):

        if read.is_secondary:
            continue

        if remove_duplicates and read.is_duplicate:
            continue

        if read.flag >= 2048:
            continue

        sequence = read.seq

        try:
            read_length = len(sequence)
        except TypeError:
            continue

        n_count = sequence.count('N')

        if float(n_count) / float(read_length) <= 0.2:
            gc_content = int(round(float(sequence.count('C') + sequence.count('G')) / float(read_length - n_count) * 100))
            gc_content_list[gc_content] += 1

        tid = read.tid
        ref_name = ''
        if tid != -1:
            ref_name = references[tid]

        if read.is_unmapped or ref_name not in chromosome_list_with_prefix or read.mapq < mapq_threshold:
            spectrumTemp2 = spectrumUnmapped
        else:
            chromosome = ref_name[chr_offset:]

            if chromosome != lastChromosome:
                chromosomeLsEnd = bands_list[chromosome]["end"]
                chromosomeLsBand = bands_list[chromosome]["band_name"]
                lastChromosome = chromosome
                i = 0
                band = chromosomeLsBand[i]
                spectrumTemp = spectrum_list[chromosome][band]

            read_start_pos = read.pos

            while i < len(chromosomeLsEnd) - 1 and read_start_pos > chromosomeLsEnd[i]:
                i += 1
                band = chromosomeLsBand[i]
                spectrumTemp = spectrum_list[chromosome][band]

            spectrumTemp2 = spectrumTemp

        spectrumTemp2["reads"] += 1

        if repeat_threshold_calc == 'n':
            repeat_threshold = int(round(float(read_length) * repeat_threshold_set / 100))

        if not consecutive_flag and len(re.findall(patterns_regex_forward, sequence)) < repeat_threshold and len(re.findall(patterns_regex_reverse, sequence)) < repeat_threshold:
            continue
        elif consecutive_flag and not re.search("(" + patterns_regex_forward + "){" + str(repeat_threshold) + "}", sequence) and not re.search("(" + patterns_regex_reverse + "){" + str(repeat_threshold) + "}", sequence):
            continue

        filtered_file.write(read)

    #############################
    ### write read count file ###
    #############################

    readcount_file = open(out_dir + "/" + pid + "_readcount.tsv", "w")
    readcount_file.write("chr\tband\treads\n")

    for chromosome in chromosome_list + ["unmapped"]:
        for band in bands_list[chromosome]["band_name"]:
            readcount_file.write("%s\t%s\t%i\n" % (chromosome, band, spectrum_list[chromosome][band]["reads"]))

    #############################
    ### write gc content file ###
    #############################

    gc_content_file = open(out_dir + "/" + pid + "_" + sample + "_gc_content.tsv", "w")
    gc_content_file.write("gc_content_percent\tread_count\n")

    for gc_content in range(0, 101):
        gc_content_file.write("%i\t%i\n" % (gc_content, gc_content_list[gc_content]))

    ##########################
    ### close file handles ###
    ##########################

    bamfile.close()
    filtered_file.close()
    readcount_file.close()
    gc_content_file.close()

    ############################
    ### index filtered file  ###
    ############################

    pysam.index(filtered_file_path)

    ##################################
    ### sort filtered file by name ###
    ##################################

    try:
        subprocess.check_call("samtools sort -n " + filtered_file_path + " " + out_dir + "/" + pid + "_filtered_name_sorted", shell=True, stderr=subprocess.PIPE)
    except:
        subprocess.call("samtools sort -n " + filtered_file_path + " -o " + out_dir + "/" + pid + "_filtered_name_sorted.bam", shell=True)





# get the reverse complement of a DNA Sequence
def getReverseComplement(sequence):
  
	sequence_temp = sequence.replace("A", "1")
	sequence_temp = sequence_temp.replace("C", "2")
	sequence_temp = sequence_temp.replace("G", "3")
	sequence_temp = sequence_temp.replace("T", "4")
	
	sequence_temp2 = sequence_temp.replace("1", "T")
	sequence_temp2 = sequence_temp2.replace("2", "G")
	sequence_temp2 = sequence_temp2.replace("3", "C")
	sequence_temp2 = sequence_temp2.replace("4", "A")
	
	sequence_reverse_complement = sequence_temp2[::-1]  # reverses a string
	
	return sequence_reverse_complement

