import sys, os
import subprocess as sub

[data, ref, output] = sys.argv[1:4]
print(data)
print(ref)
print(output)
genes="A B C DRB1 DPA1 DPB1 DQA1 DQB1".split()
#~/anaconda3/envs/CookHLA/bin/Rscript
command="Rscript src/hibag.R {} {} {} {}".format(data,ref,output," ".join(genes))
print(command)
sub.call(command.split())
print(sub.call(command.split()))
