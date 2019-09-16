#!/usr/bin/env python
# coding=utf-8
__author__ = "Peng Li"
__license__ = "MIT"
__version__ = "2.0"
__email__ = "peng-li@outlook.com"
__date__ = "2019/01/17"

'''
prepare data for pipeline analysis.
'''
import sys
import os
import os.path
import argparse
import json
import csv
from glob import glob


def safe_make_dir(path):
        try:
            os.makedirs(path)
        except OSError:
            if os.path.exists(path):
			             pass
            else:
			             raise

def file_accessible(filepath, mode):
    ''' Check if a file exists and is accessible. '''
    try:
        f = open(filepath, mode)
        f.close()
    except IOError as e:
        return False

    return True

def link_fq(new_file,orig_file):
    cmd = "ln -s " + orig_file + " " + new_file
    os.system(cmd)

def write_readgroup_file(outfile,readgroup_name,sample_name,library_name,platform_unit,run_date,platform_name,sequencing_center):
    file = open(outfile,"w")

    file.write("readgroup_name" + "\t" + readgroup_name + "\n")
    file.write("sample_name" + "\t" + sample_name + "\n")
    file.write("library_name" + "\t" + library_name + "\n")
    file.write("platform_unit" + "\t" + platform_unit + "\n")
    file.write("run_date" + "\t" + run_date + "\n")
    file.write("platform_name" + "\t" + platform_name + "\n")
    file.write("sequencing_center" + "\t" + sequencing_center + "\n")

    file.close()


def process_sample_sheet(sample_sheet,seq_type_folder,out_dir):
    ''' make dir, soft link and make a tsv file for every readgroup of fastq  '''
    with open(sample_sheet, 'r', encoding='utf-8') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('readgroup'):
                pass
            else:
                lines = line.split("\t")
                sample_level_dir = out_dir + '/' + lines[3] + '/' + seq_type_folder
                safe_make_dir(sample_level_dir)
                if(seq_type_folder=="PAIREDEND"):
                    # link R1
                    new_file_R1 = sample_level_dir + '/' + lines[0] + "_R1.fastq.gz"
                    orig_file_R1 = lines[1]
                    link_fq(new_file_R1,orig_file_R1)
                    # link R2
                    new_file_R2 = sample_level_dir + '/' + lines[0] + "_R2.fastq.gz"
                    orig_file_R2 = lines[2]
                    link_fq(new_file_R2,orig_file_R2)
                    # readgroup tsv file
                    readgroup_tsv_file = sample_level_dir + '/' + lines[0]+".tsv"
                    write_readgroup_file(readgroup_tsv_file,lines[0],lines[3],lines[4],lines[5],lines[6],lines[7],lines[8])
                else:
                    new_file_R1 = sample_level_dir + '/' + lines[0] + "_R1.fastq.gz"
                    orig_file_R1 = lines[1]
                    link_fq(new_file_R1,orig_file_R1)
                    # readgroup tsv file
                    readgroup_tsv_file =  sample_level_dir + '/' + lines[0]+".tsv"
                    write_readgroup_file(readgroup_tsv_file,lines[0],lines[3],lines[4],lines[5],lines[6],lines[7],lines[8])


def main():
    '''make command line interface'''
    parser = argparse.ArgumentParser(prog = 'prepare_data', description = 'link raw fastq files and make dirs for pipeine analysis. ')
    parser.add_argument('--sample_sheet', action = "store",required=True, help = "The sample sheet file")
    parser.add_argument('--seq_type', action = "store",required=True,default="PE",choices=["PE","SE"],help = "the read type of the fastq files, PE for paired-end, SE for single-end.")
    parser.add_argument('--out_dir', action = "store",required=True, help = 'The output  directory ')
    args = parser.parse_args()
    # parse command line argument
    sample_sheet = os.path.abspath(args.sample_sheet)
    seq_type = args.seq_type
    out_dir = os.path.abspath(args.out_dir)
    # check sample_sheet exists
    if not file_accessible(sample_sheet,'r'):
        print(sample_sheet + " does not exist or can not read!")
	# make sure out dir exists, if not make it
    safe_make_dir(out_dir)
    # Change this line to match your filenames.
    seq_type_folder = "PAIREDEND" if seq_type=="PE" else "SINGLEEND"
    process_sample_sheet(sample_sheet,seq_type_folder,out_dir)
    print("Done!")

if __name__ == '__main__':
    main()
