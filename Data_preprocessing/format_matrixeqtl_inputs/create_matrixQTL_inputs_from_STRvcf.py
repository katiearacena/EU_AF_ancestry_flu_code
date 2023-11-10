'''
Authors : Saideep Gona, Katie Aracena

This script is intended to convert a joint called genotype VCF file
into input files ready for MatrixEQTL use . It creates both the SNPs_postitions.txt and genotypes.txt file.

Modifier : Onta
Onta modified this script to process the STR vcf output from HipSTR. The modified parts are:
(1) to exclude sites with no variation among the samples.
(2) to include variants with length > 1 and multi-allelic variants in the output.
(3) to make the length of the STR as the genotype.
(4) to removed sites with MAJOR allele frequency >= 0.9 (comparable to MAF <= 0.1).
'''

import sys, os, argparse
import glob
import re

cwd = os.getcwd()
print("Python Interpreter: ", sys.executable)

parser = argparse.ArgumentParser()
args = parser.parse_args()

args.vcf_dir = "./"
args.output_dir = <MATRIX_EQTL_INPUT_DIRECTORY> ## please change the path name accordingly

def quotes(inp_s):
    return '"' + inp_s + '"'

# First collect all SNPs 

sample_vcf = args.vcf_dir + "str_calls_EU_AF_filtered.vcf"
STR_list = []
fixed = []
patch_chr = []
AF_not_pass = []
not_enough = 0

# Create data structure for genotype data
with open(sample_vcf, "r") as v:
    with open(args.output_dir + "STR_genotypes.txt", "w") as g:
        for line in v:
            if line.startswith("#CHROM"):
                p_line = line.rstrip("\n").split("\t")
                sample_id = p_line[9:]
                samplenames = ['"'+x.split("_")[1]+'"' for x in sample_id]
                g.write("\t".join(samplenames)+"\n") 

            if line.startswith("#"):
                continue

            p_line = line.rstrip("\n").split("\t")
            if p_line[4] == ".":
                fixed.append(line)
                continue
            if re.search("gl", p_line[0]):
                patch_chr.append(line)
                continue
            ref_allele = []
            ref_allele.append(p_line[3])
            alt_alleles = p_line[4].split(",")
            alleles = ref_allele + alt_alleles
                  
            geno = p_line[9:]
            genotypes = []
            genotypes_set = set()
            missing = 0
            for x in geno:
                if x.startswith("."):
                    genotypes.append("-9")
                    missing += 1
                else:
                    s = x.split(":")
                    gt = s[0].split("|")
                    
                    geno_len1 = len(alleles[int(gt[0])])
                    geno_len2 = len(alleles[int(gt[1])])
                    genotypes.append(str(geno_len1 + geno_len2))
                    genotypes_set.add(str(geno_len1 + geno_len2))

            if missing > 3:
                not_enough += 1
                continue
                
            if len(genotypes_set) == 1:
                fixed.append(line)
                continue
  
            # Remove singletons
            info = p_line[7]
            AN = re.search("AN=(\d+);", info).group(1)
            ref_AC = []
            ref_AC.append(re.search("REFAC=(\d+);", info).group(1))
            alt_ACs = re.search(";AC=([\d|,]+);", info).group(1).split(",")
            ACs = ref_AC + alt_ACs

            overrepresent = 0
            for AC in ACs:
                AF = float(AC)/float(AN)
                if AF >= 0.9:
                    overrepresent = 1
            if overrepresent == 1:
                AF_not_pass.append(line)
                continue
  
            STR = "_".join([p_line[0], p_line[1], p_line[2]])
            STR_list.append(STR)
            g.write(quotes(STR)+"\t"+"\t".join(genotypes)+"\n")
            

# Write STRs to STR_position file
    with open(args.output_dir + "STR_positions.txt", "w") as sp:
        sp.write('"str" "chr" "pos"\n')
        STR_count = 0
        for STR in STR_list:
            p_STR = STR.split("_")
            write_line = " ".join([quotes(str(STR_count)), quotes(STR), quotes('chr'+p_STR[0]), p_STR[1]]) + "\n"
            sp.write(write_line)
            STR_count += 1

print("Num of fixed STRs: ", len(fixed))
print("Num of STRs on patch chromosomes: ", len(patch_chr))
print("Num of STRs not genotyped in more than 3 samples: ", not_enough)
print("Num of STRs with overrepresenting major allele", len(AF_not_pass))
print("Num of output STRs: ", len(STR_list))
