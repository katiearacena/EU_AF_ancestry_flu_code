import sys, os, argparse
import pickle

cwd = os.getcwd()
print("Python Interpreter: ", sys.executable)

parser = argparse.ArgumentParser()
parser.add_argument("inp", help="Path to file of rasqual tests")
parser.add_argument("out", help="Path to file raw counts. The features must match the rows in this table")
args = parser.parse_args()

with open(args.inp, "r") as inp:
    with open(args.out, "w") as out:
        for line in inp:
            out.write(line[3:])