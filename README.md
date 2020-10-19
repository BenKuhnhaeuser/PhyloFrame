# PhyloFrame_Calamoideae
Code for phylogenomic analyses in Kuhnhäuser et al. 2020 "A robust phylogenomic framework for the calamoid palms"

Notes:
- Software versions are specified in the methods section of the paper
- Exemplary workflow for taxon "Calamus"

## 1. Quality control of sequence data using FastQC
```fastqc Calamus_*.fastq --extract``` produces files for visual inspection\
This was conducted both before and after trimming

## 2. Remove adapters and low-quality reads using Trimmomatic
```java -jar trimmomatic-0.38.jar PE -threads 16 -phred33 -basein Calamus_R1_001.fastq Calamus_R1_001_Tpaired.fastq Calamus_R1_001_Tunpaired.fastq Calamus_R2_001_Tpaired.fastq Calamus_R2_001_Tunpaired.fastq ILLUMINACLIP:Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:1:true LEADING:3 TRAILING:3 MAXINFO:40:0.8 MINLEN:36```
- ```ILLUMINACLIP``` is used for removing adapters specified in the file "TruSeq3-PE-2.fa"
- ```LEADING:3``` and ```TRAILING:3``` removes trimming aretacts, which are indicated by very low Phred scores
- ```MAXINFO``` balances benefits of retaining longer reads against costs of retaining bases with errors
- ```MINLEN``` sets minimum allowed length for reads

## 3. Retrieve targeted genes using Hybpiper
```python HybPiper/reads_first.py -b PhyloPalms.fasta -r Calamus_R*Tpaired.fastq --prefix $name --bwa --cov_cutoff 3```
- ```-b PhyloPalms.fasta``` specifies target file
- ```--bwa``` uses BWA to map reads to target file
- ```--cov_cutoff``` specifes minimum coverage of SPADES contigs

```python HybPiper/retrieve_sequences.py PhyloPalms.fasta . dna``` to compile one file per gene containing homologous sequences across all taxa

## 4. Sequencing statistics
Get sequence lengths for taxa in "calamuslist.txt":\
```python HybPiper/get_seq_lengths.py PhyloPalms.fasta calamuslist.txt dna > seqlengths_genes.txt```

Compute statistics on retrieved genes per taxon in "calamuslist.txt":\
```python HybPiper/hybpiper_stats.py seqlengths_genes.txt calamuslist.txt > seqlengths_genes_stats.txt```



