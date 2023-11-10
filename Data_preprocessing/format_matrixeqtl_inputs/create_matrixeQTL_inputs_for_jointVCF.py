'''
Authors : Saideep Gona, Katie Aracena

This script is intended to convert a joint called genotype VCF file
into input files ready for MatrixEQTL use . It creates both the SNPs_postitions.txt and genotypes.txt file.

'''

import sys, os, argparse
import glob

cwd = os.getcwd()
print("Python Interpreter: ", sys.executable)

parser = argparse.ArgumentParser()
parser.add_argument("vcf_file", help="Input joint VCF file")
parser.add_argument("output_dir", help="Position output file")
args = parser.parse_args()


def quotes(inp_s):
    return '"' + inp_s + '"'

# First collect all SNPs 

sample_vcf = args.vcf_file
snp_list = []

samples = []

with open(sample_vcf, "r") as v:
    for line in v:
        if line[0] == "#":
            continue
        p_line = line.rstrip("\n").split("\t")

        if len(p_line[3]) > 1:
            continue
        if len(p_line[4]) > 1:
            continue

        snp = "_".join([p_line[0], p_line[1], p_line[3], p_line[4]]) 

        snp_list.append(snp)

# Write SNPs to SNP_position file

print("Num of SNPs: ", len(snp_list))

with open(args.output_dir + "/SNP_positions.txt", "w") as sp:
    sp.write('"snp" "chr" "pos"\n')
    snp_count = 0
    for snp in snp_list:
        p_snp = snp.split("_")
        write_line = " ".join([quotes(str(snp_count)), quotes(snp), quotes("chr" + p_snp[0]), p_snp[1]]) + "\n"
        sp.write(write_line)
        snp_count += 1

# Create data structure for genotype data
with open(sample_vcf, "r") as v:
    with open(args.output_dir + "/genotypes.txt", "w") as g:
        for line in v:
            if line.startswith("#CHROM"):
                p_line = line.rstrip("\n").split("\t")
                sample_id = p_line[9:]
                samplenames = ['"'+x.split("_")[1]+'"' for x in sample_id]
                g.write("\t".join(samplenames)+"\n") 

            if line.startswith("#"):
                continue

            p_line = line.rstrip("\n").split("\t")
            geno = p_line[9:]
            genotypes = []
            for x in geno:
                if "." in x:
                    genotypes.append("-9")
                else:
                    genotypes.append(str(int(x[0])+int(x[2])))

            snp = quotes("_".join([p_line[0], p_line[1], p_line[3], p_line[4]]))

            g.write(snp+"\t"+"\t".join(genotypes)+"\n") 
           
