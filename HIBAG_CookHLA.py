#########################
'''
Last Changed :  Nov 28th, 2023.

python HIBAG_CookHLA.py -i INPUT 
                        -o OUTPUT DIR
                        -t TOOLS
                        -r REFERENCES
                        -w WEIGHTS
                        -mem MEMORY
                        -hg HUMAN GENOME VERSION
'''


import os,sys
import subprocess as sub
import argparse
from src.BASH import BASH
from src.bag2vcfh import bag2vcfh
from src.vcf2plink import vcf2plink
import numpy as np




HIBAG_fit = "python src/hibag_fit.py"
HIBAG_prefit = "python src/hibag.py"
CookHLA = "python CookHLA.py"
MakeGeneticMap = "python -m MakeGeneticMap"
PostMerger = "python src/PostMerger_CookBAG.py"


def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected for fit parameter')



def HIBAG_RUN(_in, _ref, _out, _log, _fit):
    if _fit==False:
    	command = '{} {} {} {}'. format(HIBAG_prefit, _in, _ref, _out)
    elif _fit==True:
        command = '{} -i {} -r {} -o {}'.format(HIBAG_fit, _in, _ref, _out)
    print(command)
    BASH(command, _log)    
    BASH(command, _out+"CookHLA_OUT.log") 



def CookHLA_RUN(_in, _out, _ref, _hg, _mem,_bgl):
    command = '{} -i {} -ref {} -o {} -hg {}' .format(MakeGeneticMap, _in, _ref, _out+"AGM.CookHLA", _hg)
    BASH(command, _out+"AGM.CookHLA.log") 
    if _bgl=='4':
        command = '{} -i {} -ref {} -o {} -gm {} -ae {} -mem {} -hg {} -mp 9 -bgl4' \
                    .format(CookHLA, _in, _ref, _out+"CookHLA_OUT", \
                    _out+"AGM.CookHLA.mach_step.avg.clpsB", _out+"AGM.CookHLA.aver.erate", _mem, _hg)
    if _bgl=='5':
                command = '{} -i {} -ref {} -o {} -gm {} -ae {} -mem {} -hg {} -mp 9' \
                    .format(CookHLA, _in, _ref, _out+"CookHLA_OUT", \
                    _out+"AGM.CookHLA.mach_step.avg.clpsB", _out+"AGM.CookHLA.aver.erate", _mem, _hg)
    print(command)
    BASH(command, _out+"CookHLA_OUT.log") 


def HIBAG_CookHLA(__input, __output,__references,__tools,__weights,__mem,__hg, __fit,__bgl):
    
    if __input.endswith('.vcf') or __input.endswith('.vcf.gz'):
        __input=vcf2plink(__input)
    __output = __output if not __output.endswith('/') else __output.rstrip('/')
    if bool(os.path.dirname(__output)): os.makedirs(os.path.dirname(__output), exist_ok=True)
    __output_file = os.path.basename(__output)
    
    N_CookHLA = __tools.count('CookHLA')
    N_HIBAG = __tools.count('HIBAG')
    
    idx_CookHLA = np.where(['CookHLA' in x for x in __tools])[0]
    idx_HIBAG = np.where(['HIBAG' in x for x in __tools])[0]
    
    INPUT_LIST = __output+"/" + "input_list"
    idx=1
    if N_HIBAG > 0:
        for i in range(N_HIBAG):
            HIBAG_suffix = 'HIBAG_'+str(i+1)
            print("{}. HIBAG {} RUN\n".format(idx, i+1),flush=True)
            HIBAG_dir = __output+'/' + HIBAG_suffix+'/'
            os.makedirs(HIBAG_dir, exist_ok=True)
            print("\tThe log will be written to {}\n".format(HIBAG_dir+"HIBAG_OUT.log"),flush=True)
            
            HIBAG_RUN(__input, __references[idx_HIBAG[i]], HIBAG_dir+"HIBAG_OUT", HIBAG_dir+"HIBAG_OUT.log", __fit)
            command='rm {}'.format(HIBAG_dir+"HIBAG_OUT_temp.log")
            BASH(command)
            bag2vcfh (HIBAG_dir+"HIBAG_OUT")
            command='echo {} {}' .format(HIBAG_dir+"HIBAG_OUT.vcfh", __weights[idx_HIBAG[i]])
            print("The output file's name with a given weight is written in "+INPUT_LIST+".")
            BASH(command, INPUT_LIST) 
            idx+=1

    if N_CookHLA > 0: 
        os.chdir('./CookHLA')
        for i in range(N_CookHLA):
            Cook_suffix='CookHLA_'+str(i+1)
            print("{}. CookHLA {} RUN\n".format(idx,i+1),flush=True)
            cook_dir = __output+'/' + Cook_suffix+'/'
            print(cook_dir)
            command='mkdir -p '+cook_dir
            BASH(command)
            print("\tThe log will be written to {}\n".format(cook_dir+"CookHLA_OUT.log"),flush=True)
            CookHLA_RUN(__input, cook_dir, __references[idx_CookHLA[i]], __hg, __mem, __bgl)
            if __bgl=='4':
                command = 'echo {} {}' .format(__output+'/' + Cook_suffix + "/CookHLA_OUT.MHC.QC.exon2.3000.raw_imputation_out.vcf", __weights[idx_CookHLA[i]])
            elif __bgl=='5':
                command = 'echo {} {}' .format(__output+'/' + Cook_suffix + "/CookHLA_OUT.MHC.QC.exon2.0.5.raw_imputation_out.vcf", __weights[idx_CookHLA[i]])
            else:
                raise NameError('4 or 5 is required for beagle version.')
            BASH(command, INPUT_LIST,True) 
            print("The output file's name with a given weight is written in "+INPUT_LIST+".")
            idx+=1
        os.chdir('../')
    command='mkdir -p '+__output+'/'+'Merge'
    BASH(command)
    command='mkdir -p '+__output+'/'+'Michigan_1'
    BASH(command)
            
if __name__=="__main__":
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--input","-i",default='input/1958BC.hg19', help='\nINPUT (PLINK file)\n\n')
    parser.add_argument("--output","-o",default='output/1958BC.hg_test',help='\nOUTPUT (*.alleles)\n\n')
    parser.add_argument("--tool","-t",default='CookHLA CookHLA HIBAG',help="\nTools to merge\n\n", nargs='+')
    parser.add_argument("--reference","-r",default="Reference/CookHLA/AFR/1000G_REF.AFR.chr6.hg18.29mb-34mb.inT1DGC Reference/CookHLA/Pan-Asian/Pan-Asian_REF Reference/HIBAG_prefit/African-HLA4-hg19.RData"
    ,help='\nSuper Population\n\n', nargs='+')
    parser.add_argument("--weight","-w",default='0.3 0.3 0.2'
    ,help='\nSuper Population\n\n', nargs='+')
    parser.add_argument("--memory","-mem",default='16g',help='\nMemory\n\n')
    parser.add_argument("--human-genome","-hg",default='19',help="\nhg version\n\n")
    parser.add_argument("--hibag-fit", "-fit", default=False, nargs="?", type=str2bool)
    parser.add_argument("--bgl","-bgl", default='4')
    args=parser.parse_args()
    HIBAG_CookHLA(args.input,args.output,args.reference,args.tool, args.weight, args.memory,args.human_genome, args.hibag_fit, args.bgl)
    
