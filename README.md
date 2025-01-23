# Build HMM profiles from [MiBiG 4.0 database](https://mibig.secondarymetabolites.org/) 

## Preamble
Sometimes (meta)genomic assemblers have difficulty in assembling biosynthetic gene cluters (BGCs), so using tools like [antiSMASH](https://antismash.secondarymetabolites.org/#!/start) will not give you an accurate representation of the BGCs in your (meta)genome. For example, saxitoxin and anatoxin BGCs are hard to assemble, identify, and/or bin for some reason, also antiSMASH wont identify guanitoxin BGCs. This kinda sucks :( . One method to get around these limitations is to map the raw reads to known BGCs and/or their genes. Another method is to annotate the genes in the (meta)genome. Some of these genes are in various databases, but not all genes. Becuase of this, I pulled the BGCs and the genes from MiBiG and created HMM profiles for them. Genes needed to be >50 AA long to be considered for building an HMM profile; I chose this arbitrarily.

## Building the HMM profiles

### Download and unload MiBig 4.0 database gbk files
```
wget https://dl.secondarymetabolites.org/mibig/mibig_gbk_4.0.tar.gz
tar xvf mibig_gbk_4.0.tar.gz
```
### Make our mamba env for analyses
```
mamba create -c bioconda -c defaults -c conda-forge -c biocore -n HMM_maker hmmer pip muscle

mamba activate HMM_maker

pip install requests beautifulsoup4 biopython pandas openpyxl
```
### Run the analyses
```
mamba activate HMM_maker

python HMM_Maker_V1.py mibig_gbk_4.0 results_final
```


If you use this tool please credit the hard work from the folks at MiBiG. I am not affiliated with them
