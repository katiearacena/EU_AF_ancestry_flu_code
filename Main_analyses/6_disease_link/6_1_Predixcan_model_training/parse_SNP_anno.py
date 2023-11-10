import sys
import re
import string

if len(sys.argv) != 3:
    print("Usage: python parse_dosage_rsid.py [input] [output]")
    sys.exit (1)
    
input = open(sys.argv[1], "r")
output = open(sys.argv[2], "w")


for line in input:
    if line.startswith("chromosome"):
        output.write(line)
    else:
        fields = line.split("\t")
        varID = fields[2]
        fields[6] = varID
        newline = "\t".join(fields)
        output.write(newline)
