import sys,os,glob

flu_sbatch = "Flu_footprint.sbatch"
ni_sbatch = "NI_footprint.sbatch"

full_bed = "cts_chroms.bed"

os.system("mkdir chrom_files")
os.system("rm chrom_files/*")

awk = '''awk '{print $0 >> "''' +"chrom_files/"+'''"$1".bed"}' '''+full_bed
os.system(awk)

chrom_files = glob.glob("chrom_files/*.bed")
print(chrom_files)
os.system("mkdir footprints")

for cf in chrom_files:
    prefix = cf.split("/")[1].split(".")[0]
    os.system("sbatch Flu_footprint.sbatch " + "Flu_"+prefix + " " + cf)
    os.system("sbatch NI_footprint.sbatch " + "NI_"+prefix + " " + cf)