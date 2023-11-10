import sys
import re
import string

if len(sys.argv) != 3:
    print("Usage: python parse_gene_anno.py [input] [output]")
    sys.exit (1)
    
input = open(sys.argv[1], "r")
output = open(sys.argv[2], "w")


for line in input:
    line = re.sub('"', '', line)
    if line.startswith("Gene_ID"):
        header = "\t".join(("chr", "gene_id", "gene_name", "start", "end", "gene_type")) + "\n"
        output.write(header)
    else:
        fields = line.split("\t")
        chr = re.sub("chr", "", fields[2])
        gene_id = fields[1]
        gene_name = fields[1]
        start = fields[3]
        end = fields[4].strip()
        gene_type = "protein_coding"
        newline = "\t".join((chr, gene_id, gene_name, start, end, gene_type)) + "\n"
        output.write(newline)
