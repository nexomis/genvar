## How to 

1. Clone repository

2. DLL reference

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/375/GCF_002815375.1_ASM281537v1/GCF_002815375.1_ASM281537v1_genomic.fna.gz
mv GCF_002815375.1_ASM281537v1_genomic.fna.gz ref.fa.gz
gunzip ref.fa.gz
```

3. run for each batch

```
cd batch1/
python3 ../genvar.py ../ref.fa config.yml 
cd ../
cd batch2/
python3 ../genvar.py ../ref.fa config.yml 
```
