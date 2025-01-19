# Build HMM profiles from MiBig 4.0 database

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
