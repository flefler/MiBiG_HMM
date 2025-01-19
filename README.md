# Build HMM profiles from [MiBiG 4.0 database](https://mibig.secondarymetabolites.org/) 

## Download and unload MiBig 4.0 database gbk files
```
wget https://dl.secondarymetabolites.org/mibig/mibig_gbk_4.0.tar.gz
tar xvf mibig_gbk_4.0.tar.gz
```
## Make our mamba env for analyses
```
mamba create -c bioconda -c defaults -c conda-forge -c biocore -n HMM_maker hmmer pip muscle

mamba activate HMM_maker

pip install requests beautifulsoup4 biopython pandas openpyxl
```


If you use this tool please credit the hard work from the folks at MiBiG. I am not affiliated with them
