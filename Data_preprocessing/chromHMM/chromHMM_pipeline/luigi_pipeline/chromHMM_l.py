'''
Author : Saideep Gona

This PYTHON3 script is intended for running chromHMM

PREREQS

module load python
pip install luigi

STEPS

1.) Set up your directories as follows:

output_directory/ <- Specified in config file. All pipeline outputs will go here
        |
        data/     <- Will automatically check here
            | 
            experiment_id1/      
                |
                cellmarkfiletable.txt
                alignment_file1.bam
                alignment_file2.bam
                ...
                alignment_file3.bam
            |
            experiment_id2/
                |
                cellmarkfiletable.txt
                alignment_file1.bam
                alignment_file2.bam
                ...
                alignment_file3.bam

2.) Put all associated bam files to be analyzed together in a common 
    experiment directory (e.g. "experiment_id1" above)

3.) Create a CELLMARKFILETABLE for that experiment, call it cellmarkfiletable.txt,
    and place it in the appropriate experiment directory
    (See https://www.nature.com/articles/nprot.2017.124, "Binarization")
    The goal is to specify which samples/assay each bam file pertains to

4.) Create a "luigi.cfg" config file (template is provided) and fill in the
    fields. Place it in a reasonable location (preferably close to the pipeline
    output directories)

5.) Copy the luigi_wrapper.sh file and place it in the same directory next to 
    the "luigi.cfg" file

6.) To run:
    bash luigi_wrapper RunChromHMM --experiment-id *Experiment Name*

    
'''

import sys, os, argparse
from glob import glob
from datetime import datetime
from os.path import join

import linecache

import luigi
from luigi.contrib.external_program import ExternalProgramTask
from luigi.util import requires
from bioluigi.scheduled_external_program import ScheduledExternalProgramTask
import shlex

import threading
import subprocess
import time

# Boilerplate
cwd = os.getcwd()           # Store current working directory
print("Current Working Directory: ", cwd)
dtime = str(datetime.now())
print("Time of Running ", dtime)
print("Python Interpreter: ", sys.executable)

#####################################################################################
# CONFIG AND PREP INPUT FILES
#####################################################################################

class CHROMHMM(luigi.Config):

    GENOME = luigi.Parameter()

    SINGLE_END = luigi.Parameter()

    CHROMHMM_DIR = luigi.Parameter()
    NUM_STATES = luigi.Parameter()
    LEARNMODEL = luigi.Parameter()
    TRAINED_MODEL = luigi.Parameter()
    BINSIZE = luigi.Parameter()

    USER = luigi.Parameter()

    OUTPUT_DIR = luigi.Parameter()

    sbatch_templates = luigi.Parameter()

cfg = CHROMHMM()

os.system("mkdir -p " + cfg.OUTPUT_DIR + "/logs")

class ProduceBams(luigi.Task):
    """
    Produce the BAM files for the experiment
    """
    experiment_id = luigi.Parameter()

    def output(self):
        return [luigi.LocalTarget(f) for f in sorted(glob(join(cfg.OUTPUT_DIR, 'data', self.experiment_id,'*.bam')))]

@requires(ProduceBams)
class ReheaderChroms(luigi.Task):
    """
    """
    experiment_id = luigi.Parameter()

    def run(self):
        command = [
            "bash",
            "/project2/lbarreiro/programs/chromHMM/luigi_pipeline/rehead.sh"
        ]

        for bam in self.input():
            cp = command[:]
            cp.append(bam.path)
            os.system(" ".join(cp))



#####################################################################################
# BINARIZE AND RUN CHROMHMM
#####################################################################################

@requires(ProduceBams)
class BinarizeBAM(ScheduledExternalProgramTask):

    experiment_id = luigi.Parameter()

    scheduler = 'slurm'
    cpus = 1
    memory = 10

    def program_args(self):
        args = [
            "java",
            "-mx8000M",
            "-jar",
            os.path.join(cfg.CHROMHMM_DIR,"ChromHMM.jar"),
            "BinarizeBam",
            "-gzip",
            os.path.join(cfg.CHROMHMM_DIR,"CHROMSIZES",cfg.GENOME + ".txt"),
            os.path.join(cfg.OUTPUT_DIR, "data", self.experiment_id),
            os.path.join(cfg.OUTPUT_DIR, "data", self.experiment_id, "cellmarkfiletable.txt"),
            os.path.dirname(self.output().path),
        ]
        return args
    
    def run(self):
        self.output().makedirs()
        os.system("touch " + self.output().path)
        return super(BinarizeBAM, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, "binary_intermediates", self.experiment_id, "complete.txt"))

@requires(BinarizeBAM)
class RunChromHMM(luigi.Task):

    experiment_id = luigi.Parameter()

    # scheduler = 'slurm'
    # cpus = 20
    # memory = 50

    def run(self):
        self.output().makedirs()
        os.system("touch " + self.output().path)
        if cfg.LEARNMODEL == "True":
            args = [
                "unset",
                "DISPLAY",
                "&&",
                "java",
                "-mx8000M",
                "-jar",
                os.path.join(cfg.CHROMHMM_DIR,"ChromHMM.jar"),
                "LearnModel",
                "-p",
                "0",
                os.path.dirname(self.input().path),
                os.path.dirname(self.output().path),
                cfg.NUM_STATES,
                cfg.GENOME
            ]
        elif cfg.LEARNMODEL == "False":
            args = [
                "unset",
                "DISPLAY",
                "&&",
                "java",
                "-mx8000M",
                "-jar",
                os.path.join(cfg.CHROMHMM_DIR,"ChromHMM.jar"),
                "MakeSegmentation",
                cfg.TRAINED_MODEL,
                os.path.dirname(self.input().path),
                os.path.dirname(self.output().path)
            ]

        script = '''#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50GB
'''

        with open(self.output().path, "w") as o:
            o.write(script)
            o.write("#SBATCH -e "+join(cfg.OUTPUT_DIR, "logs", "ChromHMM_"+self.experiment_id + "%J.e")+"\n")
            o.write("#SBATCH -o "+join(cfg.OUTPUT_DIR, "logs", "ChromHMM_"+self.experiment_id + "%J.o")+"\n")
            o.write("\n")
            o.write(" ".join(args)) 

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, "final_outputs", self.experiment_id, "run.sbatch"))

