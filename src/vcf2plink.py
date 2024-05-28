import os,sys
from .BASH import BASH
PLINK1="plink --noweb --silent --allow-no-sex"
PLINK2="plink2 --silent --allow-no-sex"

def chr6(_file):

    command = '{} --bfile {} --chr 6 --make-bed --out {}' \
            .format(PLINK1, _file, _file+".chr6")
    BASH(command)
    return _file+".chr6"

def duplicatePos(_bim):

    bim=[]
    dup={}
    bim_d=[]

    with open(_bim,"r") as f1:
        for i in f1:
            if i.split()[3] not in dup.keys():
                bim.append(i)
                dup[i.split()[3]]=True
            else:
                bim_d.append(i.split()[1])

    with open(_bim.replace(".bim",".toexclude"),"w") as fw:
        for i in bim_d:
            fw.write("".join(i)+"\n")

def vcf2plink(_file):
    
    if ".gz" in _file :
        command = 'gunzip {}'.format(_file)
        BASH(command)
        _file=_file.replace(".gz","")

    _out=_file.replace(".vcf","")
    command = '{} --vcf {} --make-bed --out {}' \
            .format(PLINK2, _file, _out)
    BASH(command)

    _out=chr6(_out)
    
    duplicatePos(_out+".bim")
    
    command = '{} --bfile {} --exclude {} --make-bed --out {}' \
            .format(PLINK1, _out, _out+".toexclude", _out+".dup") 
    BASH(command)

    return _out+".dup"
