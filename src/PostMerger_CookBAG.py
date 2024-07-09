
HLA=["HLA_A","HLA_B","HLA_C","HLA_DRB1","HLA_DQA1","HLA_DQB1","HLA_DPA1","HLA_DPB1"]
exon=["2","3","4"]
overlap=["3000","4000","5000"]
overlap_bgl5=["0.5","1","1.5"]

def HLA_extract_CookHLA(vcf_input):
    for i in HLA:
        os.system("grep '#CHROM' %s | sed -e 's/^#//' > %s.%s"%(vcf_input,vcf_input,i))
        os.system("grep '%s' %s >> %s.%s"%(i,vcf_input,vcf_input,i))
def HLA_extract_HIBAG(vcfh_input):
    for i in HLA:
        os.system("head -1 %s > %s.%s"%(vcfh_input,vcfh_input,i))
        os.system("grep '%s' %s >> %s.%s"%(i,vcfh_input,vcfh_input,i))
def VcfWeight(__in):
    vcf_name,gene=__in[0].split(),__in[1]
    weight=float(__in[2])

    for i in range(len(vcf_name)):
        vcf_name[i]+="."+gene

    vcf=[[] for _ in range(9)]

    for i in range(len(vcf_name)):
        with open(vcf_name[i],"r") as f1:        
            for j in f1:
                vcf[i].append(j.split())

    row=len(vcf[0])
    col=len(vcf[0][0])
    new_HLA_vcf=[[0 for i in range(col)] for j in range(row)]
    new_HLA_vcf[0][9:]=vcf[0][0][9:]
    for i in range(row): 
        new_HLA_vcf[i][:9]=vcf[0][i][:9]

    for i in range(1,row):
        for j in range(9,col):
            HLApb1,HLApb2,HLApb3=0,0,0
            for k in range(0,3):
                if HLApb1<float(vcf[k][i][j].split(":")[2].split(",")[0])+float(vcf[k][i][j].split(":")[2].split(",")[1])/2:
                    HLApb1=float(vcf[k][i][j].split(":")[2].split(",")[0])+float(vcf[k][i][j].split(":")[2].split(",")[1])/2

            for k in range(3,6):
                if HLApb2<float(vcf[k][i][j].split(":")[2].split(",")[0])+float(vcf[k][i][j].split(":")[2].split(",")[1])/2:
                    HLApb2=float(vcf[k][i][j].split(":")[2].split(",")[0])+float(vcf[k][i][j].split(":")[2].split(",")[1])/2
            
            if len(vcf[6])>i:
                for k in range(6,9):
                    if HLApb3<float(vcf[k][i][j].split(":")[2].split(",")[0])+float(vcf[k][i][j].split(":")[2].split(",")[1])/2:
                        HLApb3=float(vcf[k][i][j].split(":")[2].split(",")[0])+float(vcf[k][i][j].split(":")[2].split(",")[1])/2        
                new_HLA_vcf[i][j]=str(round((HLApb1+HLApb2+HLApb3)/3*weight,4))

            else: new_HLA_vcf[i][j]=str(round((HLApb1+HLApb2)/2*weight,4))

    with open(vcf_name[0]+"."+str(weight),"w") as f1:
        for i in new_HLA_vcf:
            f1.write("\t".join(i)+"\n")

def VcfhWeight(__in):

    vcfh_name,gene=__in[0:2]
    weight=float(__in[2])
    with open(vcfh_name+"."+gene,"r") as f1, open(vcfh_name+"."+gene+"."+str(weight),"w") as fw:
        fw.write(next(f1))
        for i in f1:
            p=i.split()
            for j in range(9,len(p)):
                p[j]=str(round(float(p[j])*weight,4))
            fw.write("\t".join(p)+"\n")

def VcfMerge(HLA_allele):
    alleles_dict={}
    alleles_pp_dict={}
    alleles_merged=[]
    for i in range(len(vcf_file)):
        alleles_dict[i]={}
        with open(vcf_file[i][0]+".HLA_"+HLA_allele+"."+vcf_file[i][1],"r") as vcf:
            FID=next(vcf).split()[9:]
            for j in FID: alleles_dict[i][j]={}
            for j in vcf:
                pp=j.split()[9:]
                HLA_vcf=j.split()[2].split("_exon")[0]
                for l in range(len(pp)):
                    alleles_dict[i][FID[l]][HLA_vcf]=pp[l]
        alleles_merged=alleles_merged+list(alleles_dict[i][FID[0]].keys())
    alleles_merged=list(sorted(set(alleles_merged)))
    if not len(alleles_merged):
        with open(output+"."+HLA_allele+".alleles","w") as fw:
            fw.write("")
        return
    ### alleles_dict[vcf_files][FID][HLA]
    alleles_answer=[[0 for i in range(8)] for j in FID]
    cnt=0
    for i in FID:
        alleles_pp_dict[i]={}
        for j in alleles_merged:
            o,pp=0,0
            for k in range(len(vcf_file)):
                 if j in alleles_dict[k][i].keys():
                        o=o+float(vcf_file[k][1])
                        pp=pp+float(alleles_dict[k][i][j])
            if o!=0:
                alleles_pp_dict[i][j]=round(pp/o,4)
        alleles_answer[cnt][0:3]=i,i,HLA_allele
        alleles_answer[cnt][5]=max(alleles_pp_dict[i].values())
        allele_1=[k for k,v in alleles_pp_dict[i].items() if v==alleles_answer[cnt][5]][0]
        mmax,ta=0,True
        for k,v in alleles_pp_dict[i].items():
            if k!=allele_1 and v>=alleles_answer[cnt][5]/2 and mmax <= v:
                mmax=v
                allele_2=k
                alleles_answer[cnt][6]=v
                alleles_answer[cnt][7]=round(v+alleles_answer[cnt][5],4)
                ta=False
        if ta:
            allele_2=allele_1
            alleles_answer[cnt][6]=alleles_answer[cnt][5]
            alleles_answer[cnt][7]=alleles_answer[cnt][5]
        alleles_answer[cnt][3]=allele_1.split("_")[-1][0:2]+","+allele_2.split("_")[-1][0:2]
        pp,qq=allele_1.split("_")[-1],allele_2.split("_")[-1]
        if HLA_allele=='DRB1':
            if allele_1.split("_")[-1]=='1454' : pp='1401'
            if allele_2.split("_")[-1]=='1454' : qq='1401'
        if alleles_answer[cnt][7] > 1:  alleles_answer[cnt][7]=1.0
        alleles_answer[cnt][4]=pp+","+qq
        cnt+=1
    with open(output+"."+HLA_allele+".alleles","w") as fw:
        for i in alleles_answer:
            fw.write(" ".join(map(str,i))+"\n")
def clear():
    
    for i in HLA: os.system("rm %s"%(output+"."+i.split("HLA_")[1]+".alleles"+"\t"))

    for i in HLA:
        for j in vcf_files.keys():
            for l in j.split():
                os.system("rm %s"%(l+"."+i+"*"))
        for j in vcfh_file:
            os.system("rm %s"%(j[0]+"."+i+"*"))

if __name__=="__main__":
     
    import os,sys,re,subprocess
    import argparse

    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--input","-i", help='input list(vcf,vcfh)',metavar='input list')
    parser.add_argument("--output","-o", help='output',metavar='Output')
    args=parser.parse_args()
 
    vcf_list,output=args.input,args.output
    vcf_file,vcfh_file=[],[]
    vcf_files={}

    if '/' in output : outdir=re.search('[\S/]*/',output).group()
    else : outdir='./'
    
    '''
    vcf_file
        
            file1, weight1
            file2, weight2
                ,,,

    vcf_files
            
            [list of file1]=weight1
            [list of file2]=weight2
                ,,,
    '''
    
    with open(vcf_list,"r") as f1:
        for i in f1:
            if 'vcfh' in i:
                vcfh_file.append(i.split()) 
            else:
                vcf_file.append(i.split())
                f=i.split()[0]
                fileinput=""
                for j in range(3):
                    for k in range(3):
                        if "2.3000" in f:
                            fileinput+=f.split("2.3000")[0]+exon[j]+"."+overlap[k]+f.split("2.3000")[1]+" "
                        else :
                            fileinput+=f.split("2.0.5")[0]+exon[j]+"."+overlap_bgl5[k]+f.split("2.0.5")[1]+" "
                vcf_files[fileinput]=i.split()[1]
    if vcf_file:
        for k in vcf_files.keys():
            for i in k.split():
                HLA_extract_CookHLA(i)
        for i in HLA:
            for k,v in vcf_files.items():
                VcfWeight([k,i,v])
    
    if vcfh_file:
        for i in vcfh_file:
            HLA_extract_HIBAG(i[0])
        for i in HLA:
            for j in vcfh_file:
                VcfhWeight([j[0],i,j[1]])
        for i in vcfh_file: vcf_file.append(i)

    for i in HLA:
        VcfMerge(i.split("HLA_")[1])
    
    alleles_list=""
    for i in HLA:
        alleles_list=alleles_list+output+"."+i.split("HLA_")[1]+".alleles"+"\t"

    os.system("cat %s > %s"%(alleles_list,output+".alleles"))
    
    clear() 
