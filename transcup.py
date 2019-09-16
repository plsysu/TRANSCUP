#!/usr/bin/env python
# coding=utf-8
__author__ = "Peng Li"
__license__ = "MIT"
__version__ = "1.0"
__email__ = "peng-li@outlook.com"
__date__ = "2019/7/27"

'''
Make a snakefile for TRANSCUP analysis.
'''

import sys
import os
import os.path
import argparse
import json
from glob import glob
import logging
import multiprocessing




CWD = os.path.dirname(os.path.realpath(__file__))

print("Current transcup.py directory: " + CWD)
print("\n")
main_pipeline_path = CWD + "/pipeline/rna_snake_v3.1.py"
default_config = CWD + "/config/database.software.resource.json"
default_config_cluster =  CWD + "/config/slurm.cluster.resource.json"



def safe_make_dir(path):
        try:
            os.makedirs(path)
        except OSError:
            if os.path.exists(path):
                         pass
            else:
                         raise

def create_logging(log_name):
        #  create a logger
        logger = logging.getLogger(log_name)
        logger.setLevel(logging.DEBUG)
        # create a handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        # define handler output format
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)

        logger.addHandler(ch)
        return logger

def file_accessible(filepath, mode):
    ''' Check if a file exists and is accessible. '''
    try:
        f = open(filepath, mode)
        f.close()
    except IOError as e:
        return False

    return True

def flat_list(list_of_lists):
    flattened = []
    for sublist in list_of_lists:
            for val in sublist:
                flattened.append(val)
    return flattened




def main():
    '''make command line interface'''
    parser = argparse.ArgumentParser(prog = 'transcup.py', description = 'Make a snakefile for TRANSCUP analysis. ')
    parser.add_argument('--fq_dir', action = "store",required=True, help = 'The parepared fastq directory by prepare_data.py ')
    parser.add_argument('--config', action = "store",required=False, default=default_config,help = 'The database/software resource config json file')
    parser.add_argument('--cluster_config', action = "store",required=False, default=default_config_cluster,help = 'The slurm cluster resource config json file')
    parser.add_argument('--out_dir', action = "store",required=True, help = 'The output directory')
    args = parser.parse_args()
    # parse command line argument
    fq_dir =  os.path.abspath(args.fq_dir)
    fq_dir = fq_dir + '/'

    config_file = os.path.abspath(args.config)
    config_file_cluster = os.path.abspath(args.cluster_config)

    if not file_accessible(config_file,'r'):
        logger.warning(config_file + " does not exist or can not read!")
        sys.exit()

    if not file_accessible(config_file_cluster,'r'):
        logger.warning(config_file_cluster + " does not exist or can not read!")
        sys.exit()

    project_dir = os.path.abspath(args.out_dir)



    logger = create_logging('TRANSCUP logger')

    safe_make_dir(project_dir)
    logger.info(project_dir + " directory was made sucessfully!")
    # output files dirs and names
    out_snakefile =  project_dir + "/Snakefile"
    out_run_cluster_sh =  project_dir + "/RUN_slurm_cluster.sh"
    out_run_local_sh =  project_dir + "/RUN_local.sh"
    out_run_SGE_cluster_sh =  project_dir + "/RUN_SGE_cluster.sh"

    # output run sh
    cluster_cmd = "sbatch --job-name {cluster.job-name} --ntasks {cluster.ntasks}  --cpus-per-task {threads} --mem {cluster.mem} --partition {cluster.partition}  --time {cluster.time} --error {cluster.error} --output {cluster.output}"
    with open(out_run_cluster_sh, 'w') as OUT_SH:
        OUT_SH.write("cluster_log=" + "./logs" + "\n")
        OUT_SH.write("mkdir -p " + "${cluster_log}" + "\n")
        OUT_SH.write("snakemake \\" + "\n")
        OUT_SH.write("--snakefile " + out_snakefile + " \\" + "\n")
        OUT_SH.write("--jobs 999 -pr --keep-going --local-cores 4 --latency-wait 120" + " \\" + "\n")
        OUT_SH.write("--configfile " +config_file + " \\" + "\n")
        OUT_SH.write("--cluster-config " +config_file_cluster + " \\" + "\n")
        OUT_SH.write("--directory " +project_dir + " \\" + "\n")
        OUT_SH.write("--cluster " + "\"" + cluster_cmd + "\"" + "\n")
        OUT_SH.write("\n")
    OUT_SH.close()
    logger.info(out_run_cluster_sh + " was created sucessfully!")

    # output run sh

    cluster_cmd_SGE = "qsub -cwd -l vf={cluster.mem} -pe smp {cluster.threads} -q all.q -V -o {cluster.output} -e {cluster.error}"
    with open(out_run_SGE_cluster_sh, 'w') as OUT_SH_SGE:
        OUT_SH_SGE.write("cluster_log=" + "./logs" + "\n")
        OUT_SH_SGE.write("mkdir -p " + "${cluster_log}" + "\n")
        OUT_SH_SGE.write("snakemake \\" + "\n")
        OUT_SH_SGE.write("--snakefile " + out_snakefile + " \\" + "\n")
        OUT_SH_SGE.write("--jobs 999 -pr --keep-going --local-cores 4 --latency-wait 120" + " \\" + "\n")
        OUT_SH_SGE.write("--configfile " +config_file + " \\" + "\n")
        OUT_SH_SGE.write("--cluster-config " +config_file_cluster + " \\" + "\n")
        OUT_SH_SGE.write("--directory " +project_dir + " \\" + "\n")
        OUT_SH_SGE.write("--cluster " + "\"" + cluster_cmd_SGE + "\"" + "\n")
        OUT_SH_SGE.write("\n")
    OUT_SH.close()
    logger.info(out_run_SGE_cluster_sh + " was created sucessfully!")


    CORES = multiprocessing.cpu_count()
    CORES_USED = str(CORES-2)
    with open(out_run_local_sh, 'w') as OUT_SH2:
        OUT_SH2.write("snakemake \\" + "\n")
        OUT_SH2.write("--snakefile " + out_snakefile + " \\" + "\n")
        OUT_SH2.write("--configfile " +config_file + " \\" + "\n")
        OUT_SH2.write("--cores " + CORES_USED + " --rerun-incomplete \n")
    OUT_SH2.close()
    logger.info(out_run_local_sh + " was created sucessfully!")

    project_dir = project_dir + '/'
    with open(out_snakefile, 'w') as OUTPROFILE:
        OUTPROFILE.write("import os, glob, sys\n")
        OUTPROFILE.write("\n")
        OUTPROFILE.write("ROOTFOLDER =  "+ "\"" + project_dir + "\"" + "\n")
        OUTPROFILE.write("FASTQFOLDER =  "+ "\"" + fq_dir + "\"" + "\n")
        OUTPROFILE.write("include:"+ "\"" + main_pipeline_path +  "\"" + "\n")
        lines = [
"""
rule run:
    input:
        ROOTFOLDER + 'complete.txt'
\n""",

                ]

        OUTPROFILE.writelines(lines)
        OUTPROFILE.close()
    logger.info(out_snakefile + " was created sucessfully!")
    logger.info("Done!")
if __name__ == '__main__':
    main()
