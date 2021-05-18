
import sys
import os
from pathlib import Path
import cassiopeia.ProcessingPipeline.process as process
import pandas as pd 
import subprocess
import cassiopeia
import numpy as np

import tqdm


#SAM_HEADER_PCT48 = "@HD	VN:1.3\n@SQ	SN:LV641_amplicon.ref	LN:271" 

def convert_sam_to_bam(sam_input, bam_output):
    """
    Converts a SAM file to BAM file.
    :param sam_input:
        Sam input file 
    :param bam_output:
        File path to write the bam output.
    :return:
        None.
    """

    cmd = "samtools view -S -b " + sam_input + " > " + bam_output

    os.system(cmd)

def call_indels(alignments, ref, output, context=True):
    """
    Given many alignments, we extract the indels by comparing the CIGAR strings in the alignments to the reference sequence. 
    :param alignments:
        Alignments provided in SAM or BAM format.
    :param ref:
        File path to the reference seqeunce, assumed to be a FASTA.
    :param output:
        Output file path.
    :param context:
        Include sequence context around indels.
    :return:
        None
    """

    perl_script = "/home/panosbino/Desktop/callAllelesLV641.pl"
    cmd = "perl " + str(perl_script) + " " + alignments + " " + ref + " " + output
    if context:
        cmd += " --context"

    cmd += " > _log.stdout 2> _log.stderr"
    print(cmd)
    p = subprocess.Popen(cmd, shell=True)
    pid, ecode = os.waitpid(p.pid, 0)

    bam_file = str(Path(output).with_suffix(".bam"))

    convert_sam_to_bam(output, bam_file) 

def align_sequences(ref, queries, outfile, gapopen=20, gapextend=1, ref_format="fasta", query_format="fastq", out_format="sam"):
    """
    Aligns many queries to a single reference sequence using EMBOSS water. By default, we assume the reference is in FASTA format, the queries 
    are all in FASTQ format, and that the output format will be a SAM file. The output file is automatically written.

    :param ref:
        File path to the reference sequence.
    :param queries:
        Queries, provided as a dataframe output from the `pickSeq` function. This will automatically be converted to a FASTQ file.
    :param gapopen:
        Gap open penalty.
    :param gapextend:
        Gap extension penalty.
    :param ref_format:
        Format of reference sequence.
    :param query_format:
        Format of query seqeunces.
    :param out_format:
        Output file format.
    :return:
        None. 
    """


    
    queries_fastq = str(Path(queries).with_suffix(".fastq"))
    collapseDF2Fastq(queries, queries_fastq)  

    cmd = "water -asequence " + ref + " -sformat1 " + ref_format  + " -bsequence " + queries_fastq + " -sformat2 " + query_format + " -gapopen " + str(gapopen) + " -gapextend " + str(gapextend) + " -outfile " + outfile + " -aformat3 " + out_format

    cmd = cmd.split(" ")

    subprocess.check_output(cmd)

    with open(outfile, "r+") as f:
        content = f.read()
        f.seek(0,0)
        f.write(SAM_HEADER + "\n" + content)#Changed by Panos


input_dir = "/proj/uppstore2018019/private/panos/10xData/10x94/outs"
home_dir = "/proj/uppstore2018019/private/panos/10xData/10x94/Cassiopeia_analysis"
ref_dir = "/proj/uppstore2018019/private/panos"
genome_bam = input_dir + "/possorted_genome_bam.bam"


#Added by Panos
ref_name = ref.split("/")[-1]
ref_name = ref_name.replace(".fa","")
with open(ref, "r+") as a_ref:
    lines = a_ref.readlines()
    Length = len(lines[1]) #Reference file has to be a one liner!!!!
#Added by Panos

SAM_HEADER = "@HD   VN:1.3\n@SQ SN:{}    LN:{}".format(ref_name, Length)

process.collapseUMIs(home_dir, genome_bam, force_sort=True)

process.pickSeq( home_dir + "/possorted_genome_bam_sorted.collapsed.txt", home_dir + "/possorted_genome_bam.picked.txt", home_dir, cell_umi_thresh=5, avg_reads_per_UMI_thresh=1, save_output=True)


align_sequences( ref_dir + "/LV641_amplicon.ref.fa", home_dir + "/possorted_genome_bam.picked.txt", home_dir + "/sw_aligned.sam")

call_indels(home_dir + "/sw_aligned.sam", ref_dir + "/LV641_amplicon.ref.fa", home_dir + "/umi_table.sam")

process.errorCorrectUMIs( home_dir + "/umi_table.bam", "", "ec_log.txt")

process.filter_molecule_table( home_dir + "/umi_table_sorted_ec.moleculeTable.txt", "moleculeTable.filtered.txt", home_dir, doublet_threshold=0.1)

process.changeCellBCID(home_dir + "/moleculeTable.filtered.txt", "", home_dir + "/moleculeTable.filtered.sample.txt")

process.call_lineage_groups( home_dir + "/test.moleculeTable.filtered.sample.txt", "alleleTable.txt", "lg_output")


#perl home/panosbino/Desktop/callAllelesLV641.pl /home/panosbino/Desktop/amp161/sw_aligned.sam /home/panosbino/Desktop/amp161/LV641_amplicon.ref.fa /home/panosbino/Desktop/amp161/umi_table.sam
