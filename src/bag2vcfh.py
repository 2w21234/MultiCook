import os,sys

def bag2vcfh(bag_input):

    HLA=["A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1"]
    bag_inputlist=[]
    bag_hla={}
    for i in HLA:
        f=bag_input+"_"+i+".bagout"
        bag_inputlist.append(f)
        bag_hla[f]="HLA_"+i

    output=bag_input+".vcfh"

    header=False
    for i in bag_inputlist:
        if not os.path.exists(i): continue 
        bag=[]
        with open(i,"r") as bagout:
            p=next(bagout).split()
            p.insert(0," ")
            bag.append(p)
            for j in bagout:
                bag.append(j.split())
    
        alleles=[]
        for j in range(1,len(bag)):
            first_allele="".join(bag[j][0].split('/')[0].split(":"))
            if first_allele == "NANANANA" or first_allele in alleles : break
            alleles.append(first_allele)

        bag_dict={}
        for j in range(1,len(bag[0])):
            ind=bag[0][j]
            bag_dict[ind]={}
            for k in alleles:
                bag_dict[ind][k]=0.0

            pp=[float(bag[k][j]) for k in range(1,len(bag))]
        
            max_1=max(pp)
            nmax_1=pp.index(max_1)
            allele1="".join(bag[nmax_1+1][0].split('/')[0].split(":"))
            allele2="".join(bag[nmax_1+1][0].split('/')[1].split(":"))
            if allele1=="NANANANA" : allele1=allele2
            bag_dict[ind][allele1]=0.5
            bag_dict[ind][allele2]=0.5
    
        if not header:
            with open(output,"w") as fw:
                for j in range(9): fw.write("0\t")
                for j in range(1,len(bag[0])): fw.write(bag[0][j]+"\t")
                fw.write("\n")
                header=True
        with open(output,"a") as fw:
            for j in alleles:
                fw.write("0\t0\t")
                fw.write(bag_hla[i]+"_"+j+"_exon\t")
                for k in range(6): fw.write("0\t")
                for k in bag_dict.keys(): fw.write(str(bag_dict[k][j])+"\t")
                fw.write("\n")

