# PhyloFrame_Calamoideae
Code for phylogenomic analyses in Kuhnhäuser et al. 2020 "A robust phylogenomic framework for the calamoid palms"

Notes:
- Software versions are specified in the methods section of the paper
- Exemplary workflow for taxon "Calamoid1" and gene "Gene1"

## 1. Quality control of sequence data
```fastqc Calamoid1_R*_001.fastq --extract``` produces files for visual inspection\
- Conduct before and after trimming

## 2. Remove adapters and low-quality reads
```java -jar trimmomatic-0.38.jar PE -threads 16 -phred33 -basein Calamoid1_R1_001.fastq Calamoid1_R1_001_Tpaired.fastq Calamoid1_R1_001_Tunpaired.fastq Calamoid1_R2_001_Tpaired.fastq Calamoid1_R2_001_Tunpaired.fastq ILLUMINACLIP:Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:1:true LEADING:3 TRAILING:3 MAXINFO:40:0.8 MINLEN:36```
- ```ILLUMINACLIP``` is used for removing adapters specified in the file "TruSeq3-PE-2.fa"
- ```LEADING:3``` and ```TRAILING:3``` removes trimming aretacts, which are indicated by very low Phred scores
- ```MAXINFO``` balances benefits of retaining longer reads against costs of retaining bases with errors
- ```MINLEN``` sets minimum allowed length for reads

## 3. Retrieve targeted exons
1) Retrieve exons per taxon\
```python HybPiper/reads_first.py -b PhyloPalm.fasta -r Calamoid1_R*Tpaired.fastq --prefix $name --bwa --cov_cutoff 3```
- ```-b PhyloPalms.fasta``` specifies target file
- ```--bwa``` uses BWA to map reads to target file
- ```--cov_cutoff``` specifes minimum coverage of SPADES contigs
2) Compile gene files containing homologous sequences across all taxa\
```python HybPiper/retrieve_sequences.py PhyloPalm.fasta . dna```
3) Produce exon list\
```for f in *.FNA; do (echo ${f/.FNA} >> genenames.txt); done```

## 4. Compute statistics for retrieval of exons
```python HybPiper/get_seq_lengths.py PhyloPalm.fasta namelist.txt dna > seqlengths_genes.txt```\
```python HybPiper/hybpiper_stats.py seqlengths_genes.txt namelist.txt > seqlengths_genes_stats.txt```
- "namelist.txt" contains taxon names

## 5. Retrieve supercontigs
1) Retrieve introns\
```python HybPiper/intronerate.py --prefix Calamoid1```\
```python HybPiper/retrieve_sequences.py PhyloPalm.fasta . intron```

2) Retrieve supercontigs (combining exons and introns)\
```python HybPiper/retrieve_sequences.py PhyloPalm.fasta . supercontig```

## 6. Exclude paralogs from analysis
### Identify paralogs
1) Write paralogs to individual file for the taxon\
```python HybPiper/paralog_investigator.py Calamoid1 2> Calamoid1_paralogs.txt```
2) Combine paralog files of different sequences\
```cat *paralogs.txt > paralog_summary.txt```
3) Return each occurrence of genes starting with HEY or EGU (as all genes in the PhyloPalms.fasta do)\
```grep -ow 'HEY\w*\|EGU\w*' paralog_summary.txt > paralog_temp.txt``` 
4) Count number of occurrences per paralog\
```tr -c '[:alnum:]' '[\n*]' < paralogs_temp.txt | sort | uniq -c | sort -nr > paralog_count.txt```
5) Produce final list of paralogs\
```grep -ow 'HEY\w*\|EGU\w*' paralog_count.txt > paralogs.txt``` 

### Make ortholog list (= genes - paralogs)
```grep -Fv -f paralogs.txt genenames.txt > orthologs.txt``` creates new ortholog list excluding paralogs

## 7. Multiple sequence alignment
### Align genes individually
```mafft --thread 4 --localpair --adjustdirectionaccurately --maxiterate 1000 Gene1.FNA > Gene1_aligned.fasta```
- use same parameters for supercontigs

### Trim alignments
1) Use automated1 algorithm in trimAl, which is optimized for ML phylogeny reconstruction:\
```trimal -in Gene1_aligned.fasta -out Gene1_aligned_trimmed_temp1.fasta -automated1```
2) Remove sites with >=80% gaps from automated1-trimmed alignment\
```trimal -in Gene1_aligned_trimmed_temp1.fasta -out Gene1_aligned_trimmed_temp2.fasta -gt 0.2```
3) Remove sequences covering less than 10 % of the alignment length\
```trimal -in Gene1_aligned_trimmed_temp2.fasta -out Gene1_aligned_trimmed.fasta -resoverlap 0.50 -seqoverlap 10```

### Remove uninformative orthologs 
The following alignments were excluded after visual inspection using Geneious:
- empty: HEY111
- sequences only for outgroup: HEY363, EGU105041872, EGU105055023
- only 3 ingroup sequences (+4 outgroup seq), few informative sites: HEY168

## 8. Gene trees
```raxmlHPC-PTHREADS -T 3 -m GTRGAMMA -f a -p 12345 -x 12345 -# 100 -k -s Gene1_aligned_trimmed.fasta -n Gene1.tree```
- use same parameters for supercontigs

## 9. Species trees
### Coalescence

### Concatenation

## 10. Rooting the tree

## 11. Comparative analyses
