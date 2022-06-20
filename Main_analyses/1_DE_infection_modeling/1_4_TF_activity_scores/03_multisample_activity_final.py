'''
Author : Saideep Gona

Run activity scores on each sample individually using the common 
footprints set. Then compute meta pvals from the individual results
'''

import sys, os, argparse
import glob
import scipy
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

cwd = os.getcwd()
print("Python Interpreter: ", sys.executable)

# Path containing individual-specific bam files
bams = "/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/data/ATACseq_alignment"

# Path to meta footprints file
footprints = "change_footprints.bed"

# Path to "final.sbatch" script
sbatch_t = "final.sbatch"

# Path to output directory
out_root = "outputs"

# Tag to add uniqueness to naming
tag = "change"


os.system("mkdir -p " + out_root)


remaining = [
    "AF04",
    "AF06",
    "AF10",
    "AF14",
    "EU09",
    "EU17",
    "EU29",
    "EU37",
    "EU43",
    "EU47"
]

exclude = [
    "AF10",
    "EU47",
    "EU43"
]
inds = set()

# Set list of individuals to run analysis on
all_individuals = []


# Run activity analysis on samples individually

for ind in all_individuals:

    if ind in inds:
        continue
    else:
        inds.add(ind)

    with open(sbatch_t, "r") as sb_t:
        
        outdir_cur = out_root + ind
        os.system("mkdir " + outdir_cur)
        os.system("mkdir " + os.path.join("sbatch",ind))

        sbatch_temp = sb_t.read()
        print(ind)
        with open(os.path.join("sbatch",ind,ind+"_"+tag+"_act.sbatch"),"w") as out:
            sbatch_temp = sbatch_temp.replace("change/change/change_footprints_mpbs.bed",mpbs)
            sbatch_temp = sbatch_temp.replace("$PWD/OUTPUTS",outdir_cur)
            sbatch_temp = sbatch_temp.replace("NI_BAM",os.path.join(bams,ind+"_NI.sorted.dup.bam"))
            sbatch_temp = sbatch_temp.replace("FLU_BAM",os.path.join(bams,ind+"_Flu.sorted.dup.bam"))
            sbatch_temp = sbatch_temp.replace("--output-location=$PWD/OUTPUTS/IND","--output-location="+outdir_cur)

            out.write(sbatch_temp)

        os.system("sbatch " + os.path.join("sbatch",ind,ind+"_"+tag+"_act.sbatch"))

# Collect pvalues from activity output

tfs = []
inds = set()
for ind in all_individuals:

    if ind in exclude:
        continue
    if ind in inds:
        continue
    else:
        inds.add(ind)
    outdir_cur = out_root + ind

    with open(outdir_cur+"/differential_statistics.txt", "r") as stats:
        c=0
        for line in stats:
            if c == 0:
                c+=1
                continue
            p_line = line.strip().split("\t")
            tfs.append(p_line[0])

    break

cols = 0
inds = set()
for ind in all_individuals:

    if ind in exclude:
        continue
    if ind in inds:
        continue
    else:
        inds.add(ind)

    cols += 1

all_p = np.zeros((len(tfs),cols))
all_a = np.zeros((len(tfs),cols))

t = 0
tf_inds = {}
for tf in tfs:
    tf_inds[tf] = t
    t+=1

e=0
inds = set()
with open(out_root+"/samples.txt", "w") as osamp:

    for ind in all_individuals:

        if ind in exclude:
            continue
        if ind in inds:
            continue
        else:
            inds.add(ind)
        osamp.write(ind + "\n")

        outdir_cur = out_root + ind

        with open(outdir_cur+"/differential_statistics.txt", "r") as stats:
            k=0
            for line in stats:
                if k == 0:
                    k+=1
                    continue
                p_line = line.strip().split("\t")
                tf = p_line[0]
                c = tf_inds[tf]
                all_p[c,e] = float(p_line[8])
                all_a[c,e] = p_line[6]
                c+=1

        e+=1

# Compute meta pvalues

outs = np.zeros((len(tfs),3))
with open(out_root+"/tfs.txt", "w") as otfs:
    for x in range(len(tfs)):
        # print(scipy.stats.combine_pvalues(all_p[x,:]))
        outs[x,1] = scipy.stats.combine_pvalues(all_p[x,:])[1]
        outs[x,0] = np.mean(all_a[x,:])
        outs[x,2] = np.var(all_a[x,:])
        otfs.write(tfs[x] + "\n")

print(all_a)

np.savetxt(out_root + "/multisample_act.txt", outs)
np.savetxt(out_root + "/multisample_a.txt", all_a)
np.savetxt(out_root + "/multisample_p.txt", all_p)
