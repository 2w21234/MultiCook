
# MultiCook
## (1) Introduction

**Human Leukocyte Antigen (HLA)** molecules are produced by genes within Major Histocompatibility Complex (MHC). The identification of HLA genotype is challenging due to the complex genetic structure of the MHC region and is often costly. Fortunately, recent computational methods, such as CookHLA and Michigan imputation server have made it possible to impute HLA genotypes using inexpensive single nucleotide polymorphism (SNP) markers. While CookHLA performs well with ethnicity-matched panels, the availability of large-sized panels specific to each ethnicity still remains limited. The Michigan HLA imputation server leverages a large multi-ethnic reference panel. However, because of the heterogeneity of HLA allele composition across ethnicities, the efficacy of utilizing a large multi-ethnic reference panel for HLA imputation needs to be investigated.   

To achieve a better HLA imputation, we need a reference group with a similar ethnicity to the sample. But we also need a reference which is large enough to contain diverse individuals or sub-ethnic groups to increase the chance that any of them would be similar to the target sample. Addressing these issues, we introduce MultiCook which utilizes multiple imputation results from diverse imputation tools and reference panels. MultiCook calculates the posterior probability for each imputation result and subsequently linearly combines these posterior probabilities.


We introduce ```MultiCook``` which utilizes multiple imputation results from diverse imputation tools (CookHLA, Michigan Imputation Server, HIBAG) and reference panels. 
MultiCook calculates the posterior probability for each imputation result and subsequently linearly combines these posterior probabilities. 


![image](https://github.com/2w21234/MultiCook/assets/37434378/69c71ace-2502-4d79-b99b-0bf12e7fdb4b) 


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
conda create --file MultiCook.yaml
conda activate MultiCook
```

## (3) Running code
MultiCook combines the results of CookHLA, Michigan server imputation and HIBAG. 
MultiCook provides a convenient script(**Imputation_single.py**) for running single-panel-based CookHLA and HIBAG simultaneously.

```
python Imputation_single.py -i input/HapMap -t CookHLA HIBAG -o output/HapMap -r Reference/CookHLA/EUR/1000G_REF.EUR.chr6.hg18.29mb-34mb.inT1DGC Reference/HIBAG_prefit/European-HLA4-hg19.RData -w 0.5 0.5 -m 16g -hg 19 
```
**Imputation_single.py** takes 7 arguments.  
-i : the path to the input (target) data in plink binary format  
-t : names of tools (CookHLA or HIBAG)  
-o : the name of folder to contain the results of the imputaion
-r : paths to the reference panels  
-w : weights for the reference panels  
-m : the memory allocated for running CookHLA  
-hg : the version of human genome assembly of input data (hg19 or hg18)  




After running **Imputation_single.py**, a file named **input_list** containing the paths to the results from single-panel-based tools with their weights is generated.
**input_list** is the input for the merger (**Merge.py**) of the results and can be manually 

```
python Merge.py -i /data01/hakin/tmp/MultiCook/output/HapMap/input_list -o /data01/hakin/tmp/MultiCook/output/HapMap/Merge
```



