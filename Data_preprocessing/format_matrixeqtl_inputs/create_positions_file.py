'''
Author : Saideep Gona

This script converts a counts table into a position file for matrixQTL
'''

import sys, os, argparse

cwd = os.getcwd()
print("Python Interpreter: ", sys.executable)

parser = argparse.ArgumentParser()
parser.add_argument("counts_table", help="Directory containing input files")
parser.add_argument("position_file", help="Position output file")
args = parser.parse_args()

def quotes(x):
    return '"' + x + '"'

header = '"Gene_ID"\t"chromosome"\t"S1"\t"S2"\n'

count = 0

with open(args.counts_table, "r") as ct:
    with open(args.position_file, "w") as out:
        out.write(header)
        for line in ct:
            if count == 0:
                count += 1
                continue

            p_line = line.rstrip("\n").split("\t")
            gene_id = p_line[0]
            chrom = "_".join(gene_id.split("_")[0:-2])
            s1 = gene_id.split("_")[-2]
            s2 = gene_id.split("_")[-1]

            new_line = [str(count), gene_id, chrom, s1, s2]
            new_line[0] = quotes(new_line[0])
            new_line[1] = quotes(new_line[1])
            new_line[2] = quotes(new_line[2])
            new_line = "\t".join(new_line) + "\n"

            out.write(new_line)

            count += 1
