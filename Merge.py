

import os, sys, re, subprocess
import argparse
import numpy as np
from pathlib import Path
import glob

HLA = ["HLA_A", "HLA_B", "HLA_C", "HLA_DRB1", "HLA_DQA1", "HLA_DQB1", "HLA_DPA1", "HLA_DPB1"]
exon = ["2", "3", "4"]
overlap = ["3000", "4000", "5000"]
overlap_bgl5 = ["0.5", "1", "1.5"]

def HLA_extract_CookHLA(vcf_input):
    for i in HLA:
        vcf_path = Path(vcf_input)
        os.system(f"grep '#CHROM' {vcf_path} | sed -e 's/^#//' > {vcf_path}.{i}")
        os.system(f"grep '{i}' {vcf_path} >> {vcf_path}.{i}")

def HLA_extract_HIBAG(vcfh_input):
    for i in HLA:
        vcfh_path = Path(vcfh_input)
        os.system(f"head -1 {vcfh_path} > {vcfh_path}.{i}")
        os.system(f"grep '{i}' {vcfh_path} >> {vcfh_path}.{i}")

def HLA_extract_Michigan(vcf_input_M):
    for i in HLA:
        vcf_path = Path(vcf_input_M)
        os.system(f"grep '#CHROM' {vcf_path} | sed -e 's/^#//' > {vcf_path}.{i}")
        os.system(f"grep '{i}' {vcf_path} >> {vcf_path}.{i}")

def VcfWeight(__in):
    vcf_name, gene = __in[0].split(), __in[1]
    weight = float(__in[2])
    for i in range(len(vcf_name)):
        vcf_name[i] = f"{vcf_name[i]}.{gene}"
    vcf = [[] for _ in range(9)]

    for i in range(len(vcf_name)):
        with open(vcf_name[i], "r") as f1:
            for j in f1:
                vcf[i].append(j.split())

    row = len(vcf[0])
    col = len(vcf[0][0])
    new_HLA_vcf = [[0 for _ in range(col)] for _ in range(row)]
    new_HLA_vcf[0][9:] = vcf[0][0][9:]
    for i in range(row):
        new_HLA_vcf[i][:9] = vcf[0][i][:9]
    for i in range(1, row):
        for j in range(9, col):
            HLApb1, HLApb2, HLApb3 = 0, 0, 0
            for k in range(0, 3):
                if HLApb1 < float(vcf[k][i][j].split(":")[2].split(",")[0]) + float(vcf[k][i][j].split(":")[2].split(",")[1]) / 2:
                    HLApb1 = float(vcf[k][i][j].split(":")[2].split(",")[0]) + float(vcf[k][i][j].split(":")[2].split(",")[1]) / 2
            for k in range(3, 6):
                if HLApb2 < float(vcf[k][i][j].split(":")[2].split(",")[0]) + float(vcf[k][i][j].split(":")[2].split(",")[1]) / 2:
                    HLApb2 = float(vcf[k][i][j].split(":")[2].split(",")[0]) + float(vcf[k][i][j].split(":")[2].split(",")[1]) / 2

            if len(vcf[6]) > i:
                for k in range(6, 9):
                    if HLApb3 < float(vcf[k][i][j].split(":")[2].split(",")[0]) + float(vcf[k][i][j].split(":")[2].split(",")[1]) / 2:
                        HLApb3 = float(vcf[k][i][j].split(":")[2].split(",")[0]) + float(vcf[k][i][j].split(":")[2].split(",")[1]) / 2
                new_HLA_vcf[i][j] = str(round((HLApb1 + HLApb2 + HLApb3) / 3 * weight, 4))
            else:
                new_HLA_vcf[i][j] = str(round((HLApb1 + HLApb2) / 2 * weight, 4))

    with open(f"{vcf_name[0]}.{weight}", "w") as f1:
        for i in new_HLA_vcf:
            f1.write("\t".join(i) + "\n")

def VcfhWeight(__in):
    vcfh_name, gene = __in[0:2]
    weight = float(__in[2])
    with open(f"{vcfh_name}.{gene}", "r") as f1, open(f"{vcfh_name}.{gene}.{weight}", "w") as fw:
        fw.write(next(f1))
        for i in f1:
            p = i.split()
            for j in range(9, len(p)):
                p[j] = str(round(float(p[j]) * weight, 4))
            fw.write("\t".join(p) + "\n")

def VcfWeight_Michigan(__in):
    vcf_name, gene = __in[0], __in[1]
    weight = float(__in[2])
    vcf_path = Path(f"{vcf_name}.{gene}")
    with open(vcf_path, "r") as f1:
        lines = f1.readlines()
    header = lines[0]
    data_lines = lines[1:]
    new_data = []
    for line in data_lines:
        parts = line.strip().split()
        for j in range(9, len(parts)):
            try:
                parts[j] = str(round(float(parts[j]) * weight, 4))
            except ValueError:
                continue  # Skip values that cannot be converted to float
        new_data.append("\t".join(parts))
    with open(f"{vcf_name}.{gene}.{weight}", "w") as f2:
        f2.write(header)
        f2.write("\n".join(new_data) + "\n")

def VcfMerge(HLA_allele, output, vcf_file):
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)  # Ensure the output directory exists

    merged_data = {}
    for vcf in vcf_file:
        with open(f"{vcf[0]}.HLA_{HLA_allele}.{vcf[1]}", "r") as f:
            header = next(f).strip()
            for line in f:
                parts = line.strip().split()
                key = tuple(parts[:9])
                values = list(map(float, parts[9:]))
                if key not in merged_data:
                    merged_data[key] = np.array(values) * float(vcf[1])
                else:
                    merged_data[key] += np.array(values) * float(vcf[1])

    with open(f"{output_path}.{HLA_allele}.alleles", "w") as f:
        f.write(header + "\n")
        for key, values in merged_data.items():
            normalized_values = values / np.sum([float(v[1]) for v in vcf_file])
            f.write("\t".join(key) + "\t" + "\t".join(map(lambda x: str(round(x, 4)), normalized_values)) + "\n")


def clear(__output):
    for i in HLA: os.system("rm %s"%(__output+"."+i.split("HLA_")[1]+".alleles"+"\t"))

    for i in HLA:
        for j in vcf_files.keys():
            for l in j.split():
                os.system("rm %s"%(l+"."+i+"*"))
        for j in vcfh_file:
            os.system("rm %s"%(j[0]+"."+i+"*"))
        for j in vcf_file_M:
            os.system("rm %s"%(j[0]+"."+i+"*"))


            
if __name__=="__main__":


    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--input","-i",help='input list(vcf,vcfh)',metavar='input list')
    parser.add_argument("--output","-o",help='output',metavar='output')
    args=parser.parse_args()
    #args = parser.parse_args(args=[])
    #output = args.output
    
    vcf_list,output=args.input,args.output
    vcf_file,vcfh_file=[],[]
    vcf_file_M=[] # Michigan
    vcf_files={}
    
    #if '/' in output : outdir=re.search('[\S/]*/',output).group()
    #else : outdir='./'
       
    with open(vcf_list,"r") as f1:
        weights=[]
        for j in f1:
            #print(j)
            weights.append(float(j.split()[1]))
        weights=np.array(weights)
        weights=weights/np.sum(weights)
        #print(weights)

        
    with open(vcf_list,"r") as f1:
        l=0
        m=0
        for i in f1:
            #print(i)
            #print(i)
            #print([i.split()[0], weights[k]])
            if 'vcfh' in i:
                #print([i.split()[0], weights[l]])
                vcfh_file.append([i.split()[0], str(weights[l])]) 
                l += 1
            elif 'exon' in i:
                vcf_file.append([i.split()[0], str(weights[l])])
                #l += 1
                f=i.split()[0]
                fileinput=""
                for j in range(3):
                    for k in range(3):
                        if "2.3000" in f:
                            fileinput+=f.split("2.3000")[0]+exon[j]+"."+overlap[k]+f.split("2.3000")[1]+" "
                        else :
                            fileinput+=f.split("2.0.5")[0]+exon[j]+"."+overlap_bgl5[k]+f.split("2.0.5")[1]+" "
                #print(i.split()[1])
                vcf_files[fileinput]=weights[l]
                l +=1
            else:
                vcf_file_M.append([i.split()[0], str(weights[l])])
                l += 1
    #print(vcf_file)
    #print(vcfh_file)
    #print(vcf_file_M)
    #print(vcf_file)
    #print(vcfh_file)
    #print(vcf_file_M)
    if vcf_file:
        for k in vcf_files.keys():
            for i in k.split():
                HLA_extract_CookHLA(i)
        for i in HLA:
            for k,v in vcf_files.items():
                #print(j,i)
                #print(j[0])
                #print(j[1])
                #print([i,k,v])
                VcfWeight([k,i,v])
                
    if vcfh_file:
        for i in vcfh_file:
            HLA_extract_HIBAG(i[0])
        for i in HLA:
            for j in vcfh_file:
                VcfhWeight([j[0],i,j[1]])
        for i in vcfh_file: vcf_file.append(i)

    if vcf_file_M:
        HLA_extract_Michigan(vcf_file_M[0][0])
        
        for i in HLA:
            for j in vcf_file_M:
                #print(j,i)
                VcfWeight_Michigan([j[0],i,j[1]])
        for i in vcf_file_M: vcf_file.append(i)
    
    #print(vcf_file)
    print(output)
    for i in HLA:
        VcfMerge(i.split("HLA_")[1], output, vcf_file)
    
    command_sub = [output+'.'+x.split('_')[1]+".alleles " for x in HLA]
    command_sub_1 = ''
    for i in command_sub:
        command_sub_1 += i+' '
    command = 'cat '+command_sub_1+ " > " +output+".all.alleles"
    #print(command)
    os.system(command)  
    clear(output)
    
    print('Finished!')
    print('The results are in '+output+'.all.alleles')
    
