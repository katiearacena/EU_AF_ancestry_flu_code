inpu = "/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/HINT/outputs/meta/cts.bed"
outpu = "/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/HINT/outputs/meta/cts_chroms.bed"

with open(inpu, "r") as inp:
    with open(outpu, "w") as out:
        for line in inp:
            p_line = line.strip().split("\t") 
            if len(p_line[0].split("_")) > 1:
                continue
            else:
                out.write(line)