import sys, os
import subprocess as sub
import argparse


parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--fit", "-fit", default=False)
parser.add_argument("--ref-list", "-r")
parser.add_argument("--input", "-i")
parser.add_argument("--output", "-o")
args=parser.parse_args()

data = args.input
output = args.output
ref_list = args.ref_list
#[data, output] = sys.argv[1:3]
# ref_list="./Reference/HIBAG_fit/REF_LIST"
print('HI')
print(ref_list)
with open(ref_list ,"r") as f1:
    for i in f1:
        gene=i.split()[0]
        ref=i.split()[1]
        command="Rscript src/hibag_fit.R {} {} {} {}".format(data,ref,output,gene)
        print(command)
        sub.call(command.split())		

