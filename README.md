
# MultiCook
## (1) Introduction

**Human Leukocyte Antigen (HLA)** molecules are produced by genes within Major Histocompatibility Complex (MHC). The identification of HLA genotype is challenging due to the complex genetic structure of the MHC region and is often costly. Fortunately, recent computational methods, such as CookHLA and Michigan imputation server have made it possible to impute HLA genotypes using inexpensive single nucleotide polymorphism (SNP) markers. While CookHLA performs well with ethnicity-matched panels, the availability of large-sized panels specific to each ethnicity still remains limited. The Michigan HLA imputation server leverages a large multi-ethnic reference panel. However, because of the heterogeneity of HLA allele composition across ethnicities, the efficacy of utilizing a large multi-ethnic reference panel for HLA imputation needs to be investigated.   

To achieve a better HLA imputation, we need a reference group with a similar ethnicity to the sample. But we also need a reference which is large enough to contain diverse individuals or sub-ethnic groups to increase the chance that any of them would be similar to the target sample. Addressing these issues, we introduce MultiCook which utilizes multiple imputation results from diverse imputation tools and reference panels. MultiCook calculates the posterior probability for each imputation result and subsequently linearly combines these posterior probabilities.


We introduce ```MultiCook``` which utilizes multiple imputation results from diverse imputation tools (CookHLA, Michigan Imputation Server, HIBAG) and reference panels. 
MultiCook calculates the posterior probability for each imputation result and subsequently linearly combines these posterior probabilities. 


![image](https://github.com/2w21234/MultiCook/assets/37434378/69c71ace-2502-4d79-b99b-0bf12e7fdb4b) 

<br/>
<br/>

## (2) Installation
After clonning the anaconda virtual environment of MultiCook, create and activate the environment.  
The yaml file of MultiCook is based on  
> python v3.9  
> plink v1.90  
> CookHLA v1.0.1  
> HIBAG v1.34.0  


```
git clone https://github.com/2w21234/MultiCook.git
cd MultiCook
conda env create --file MultiCook.yaml
conda activate MultiCook
```
<br/>
<br/>

## (3) Running code
MultiCook combines the results of CookHLA, Michigan server imputation and HIBAG. 
MultiCook provides a convenient script (**Imputation_single.py**) for running single-panel-based CookHLA and HIBAG simultaneously.
<br/>

```
python Imputation_single.py -i input/HapMap -t CookHLA HIBAG -o output/HapMap -r Reference/CookHLA/EUR/1000G_REF.EUR.chr6.hg18.29mb-34mb.inT1DGC Reference/HIBAG_prefit/European-HLA4-hg19.RData -w 0.5 0.5 -m 16g -hg 19 
```
**Imputation_single.py** takes 7 arguments.  
> -i : the path to the input (target) data in plink binary format  
> -t : names of tools (CookHLA or HIBAG)  
> -o : the name of folder to contain the results of the imputaion  
> -r : paths to the reference panels  
> -w : weights for the reference panels  
> -m : the memory allocated for running CookHLA  
> -hg : the version of human genome assembly of input data (hg19 or hg18)  
<br/>


After running **Imputation_single.py**, a file named **input_list** containing the paths to the results from single-panel-based tools with their weights is generated.  

**input_list** is the input file for the merger (**Merge.py**) of the results. 

For CookHLA, the path to 'CookHLA_OUT.MHC.QC.exon2.3000.raw_imputation_out.vcf' must be included.  

For HIBAG, the path to 'HIBAG_OUT.vcfh' must be included.  

Users who want to merge the Michigan impuation server's reult can manually edit **input_list** to add the path to the reult ('chr6.dose.vcf') from the Michigan imputation server with its weight.  

After running Michigan imputation server, a file names 'chr6.dose.vcf' which is default name of Michigan server is generated.
To merge this result within ```MultiCook```, users should insert the path of 'chr6.dose.vcf' with its weight.  
We provide the result('chr6.dose.vcf.gz') of Michigan server, which needs to be unzipped (```gunzip chr6.dose.vcf.gz```).



```
python Merge.py -i output/HapMap/input_list -o output/HapMap/Merge/result
```


**Merge.py** takes 2 arguments.  
> -i : the path to **input_list**  
> -o : the path to the output of the merger  


After running **Merge.py**, **result.all.alleles** containing the predicted HLA allele pairs is generated.  

<br/>

CookHLA provides the code to measure the imputation accuracy.  

Under the CookHLA folder (**MultiCook/CookHLA/**), given the answer file, users can measure the accuracy as following:  

```
python -m measureAcc ../input/HapMap.answer.alleles  ../output/HapMap/Merge/result.all.alleles ../output/HapMap/Merge/result
```


Michigan imputation server requires **vcf** file, which is different from **plink bfile** of **CookHLA** or **HIBAG**.  
To merge the result of Michigan imputation server, the user should convert the format.  
This conversion should be done under a directory which includes input files(**HapMap.bed, HapMap.bim, HapMap.fam**).  
The following codes are for a conversion of toy example **plink bfile** (Refer to 'https://imputationserver.readthedocs.io/en/latest/prepare-your-data/' for details).  
Under the folder (**MultiCook/input**),

```
wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
unzip HRC-1000G-check-bim-v4.2.7.zip

plink --freq --bfile HapMap --out HapMap
perl HRC-1000G-check-bim.pl -b HapMap.bim -f HapMap.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
sh Run-plink.sh
plink --bfile HapMap-updated-chr6 --recode vcf --out HapMap
bgzip -c HapMap.vcf > HapMap.vcf.gz
rm HapMap-updated*
```
Then, **HapMap.vcf.gz** is generated for the input of Michigan imputation server.

<br/>
<br/>




## (4) References
**Kim H, Lim H and Han B. MultiCook: a tool that improves accuracy of HLA imputation by combining probabilities from multiple reference panels**

Cook S, Choi W, Lim H, et al. Accurate imputation of human leukocyte antigens with CookHLA. Nature Communications. Feb 2021;12(1)1264. doi:10.1038/s41467-021-21541-5  

Luo Y, Kanai M, Choi W, et al. A high-resolution HLA reference panel capturing global population diversity enables multi-ancestry fine-mapping in HIV host response. Nat Genet. Oct 2021;53(10):1504-1516. doi:10.1038/s41588-021-00935-7  

Zheng X, Shen J, Cox C, et al. HIBAG--HLA genotype imputation with attribute bagging. Pharmacogenomics J. Apr 2014;14(2):192-200. doi:10.1038/tpj.2013.18  

<br/>
<br/>

## (5) License
The MultiCook Software is freely available for non-commercial academic research use. For other usage, one must contact Buhm Han (BH) at buhm.han@snu.ac.kr (patent pending). WE (Hakin Kim, Hyunjoon Lim and BH) MAKE NO REPRESENTATIONS OR WARRANTIES WHATSOEVER, EITHER EXPRESS OR IMPLIED, WITH RESPECT TO THE CODE PROVIDED HERE UNDER. IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE WITH RESPECT TO CODE ARE EXPRESSLY DISCLAIMED. THE CODE IS FURNISHED "AS IS" AND "WITH ALL FAULTS" AND DOWNLOADING OR USING THE CODE IS UNDERTAKEN AT YOUR OWN RISK. TO THE FULLEST EXTENT ALLOWED BY APPLICABLE LAW, IN NO EVENT SHALL WE BE LIABLE, WHETHER IN CONTRACT, TORT, WARRANTY, OR UNDER ANY STATUTE OR ON ANY OTHER BASIS FOR SPECIAL, INCIDENTAL, INDIRECT, PUNITIVE, MULTIPLE OR CONSEQUENTIAL DAMAGES SUSTAINED BY YOU OR ANY OTHER PERSON OR ENTITY ON ACCOUNT OF USE OR POSSESSION OF THE CODE, WHETHER OR NOT FORESEEABLE AND WHETHER OR NOT WE HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES, INCLUDING WITHOUT LIMITATION DAMAGES ARISING FROM OR RELATED TO LOSS OF USE, LOSS OF DATA, DOWNTIME, OR FOR LOSS OF REVENUE, PROFITS, GOODWILL, BUSINESS OR OTHER FINANCIAL LOSS.
