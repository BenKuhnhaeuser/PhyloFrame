# PhyloFrame_Calamoideae
Scripts for phylogenomic analyses in Kuhnhäuser et al. (2021), A robust phylogenomic framework for the calamoid palms, Molecular Phylogenetics and Evolution.

Notes:
- Software versions are specified in the methods section of the paper
- Exemplary workflow for taxon "Calamoid1" and exon "Gene1"
- Identical parameters used for exons and supercontigs

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

## 9. Species trees
### Coalescence
1) Concatenate trees\
```cat RAxML_bipartitions.Gene*.tree > allTrees.tree```
2) Remove outlier branches\
```run_treeshrink.py -t allTrees.tree -q 0.05```
3) Collapse branches with bootstrap support below 10%\
```nw_ed allTrees_ts.tree 'i & b<=10' o > allTrees_ts_bs10.tree```
4) Build coalescence species tree\
```java -jar astral.5.6.3.jar -i allTrees_ts_bs10.tree -o astral.tree```
5) Root the tree
```pxrr -t astral.tree -g Nypa-fructicans-MSL30-S32,Kerriodoxa-elegans-MSL76,Asterogyne-martiana-SBL226,Ceroxylon-quindiuense-MSL17 -s > astral_rooted.tree```

### Concatenation
1) Concatenate individual gene alignments\
```AMAS concat -i Gene*_aligned_trimmed.fasta -f fasta -d dna -c 1```
2) Build concatenated species tree\
```raxmlHPC-PTHREADS -T 3 -m GTRGAMMA -f a -p 12345 -x 12345 -# autoMRE -k -s concatenated.out -n raxml_concatenated.tree```
3) Root the tree\
```pxrr -t raxml_concatenated.tree -g Nypa-fructicans-MSL30-S32,Kerriodoxa-elegans-MSL76,Asterogyne-martiana-SBL226,Ceroxylon-quindiuense-MSL17 -s > raxml_concatenated_rooted.tree```

## 10. Comparative analyses
### Tree comparisons
Load R libraries:\
```library(ape)```\
```library(phangorn)```
Read in trees:\
```concat_genes <- read.tree("concat_genes.tree")```\
```concat_super <- read.tree("concat_supercontigs.tree")```\
```astral_genes <- read.tree("astral_genes.tree")```\
```astral_super <- read.tree("astral_supercontigs.tree")```

#### Strict consensus of species trees
```consensus_all <- consensus(concat_genes, concat_super, astral_genes, astral_super)```\
```plot.phylo(consensus_all, cex = 0.7)```\
```write.tree(consensus_all, file = "consensus_all.tree")```

#### Normalized Robinson Foulds distances
```RF.dist(concat_genes, concat_super, check.labels = TRUE, normalize = TRUE, rooted = TRUE)```\
```RF.dist(concat_genes, astral_genes, check.labels = TRUE, normalize = TRUE, rooted = TRUE)```\
```RF.dist(concat_genes, astral_super, check.labels = TRUE, normalize = TRUE, rooted = TRUE)```\
```RF.dist(concat_super, astral_genes, check.labels = TRUE, normalize = TRUE, rooted = TRUE)```\
```RF.dist(concat_super, astral_super, check.labels = TRUE, normalize = TRUE, rooted = TRUE)```\
```RF.dist(astral_genes, astral_super, check.labels = TRUE, normalize = TRUE, rooted = TRUE)```

### Gene tree conflict
Exemplary command for Calameae:\
```docker run -v DiscoVista-master/:/data esayyari/discovista discoVista.py -p calamoideae/ -m 5 -a calamoideae/annotation_calameae.txt -o calamoideae/results/calameae -g Outgroup```

## 11. Statistics
### Targeted sequencing
```python HybPiper/get_seq_lengths.py PhyloPalm.fasta namelist.txt dna > seqlengths_genes.txt```\
```python HybPiper/hybpiper_stats.py seqlengths_genes.txt namelist.txt > seqlengths_genes_stats.txt```
- "namelist.txt" contains taxon names

### Alignments
```AMAS summary -f fasta -d dna -i alignment_genes_concatenated.out -o summary_alignment_genes_concatenated.txt```
- For ingroup only: first remove outgroup from concatenated alignment\
```AMAS remove -x Asterogyne-martiana-SBL226 Kerriodoxa-elegans-MSL76 Nypa-fructicans-MSL30-S32 Ceroxylon-quindiuense-MSL17 -f fasta -d dna -i alignment_genes_concatenated.out```
- For outgroup only: first remove ingroup\
```AMAS remove -x Mauritia-carana-BKL072 Raphia-hookeri-BKL088 Raphia-monbuttorum-BKL090 Raphia-textilis-BKL085 Calamus-rheedei-BKL065 Plectocomiopsis-sp-nov-discolor-RBL178 Mauritiella-aculeata-BKL010 Calamus-calospathus-BKL011 Lepidocaryum-tenue-MSL75-S14 Oncocalamus-mannii-BKL080 Myrialepis-paradoxa-BKL015 Calamus-discolor-RBL122 Raphia-farinifera-MSL20-S18 Metroxylon-sagu-MSL36-S19 Eleiodoxa-conferta-MSL73-S5 Plectocomia-elongata-RBL210 Salacca-lophospatha-BKL250 Raphia-regalis-BKL061 Calamus-radiatus-BKL176 Plectocomiopsis-geminiflora-MSL74 Calamus-peregrinus-BKL057 Plectocomia-himalayana-BKL078 Mauritiella-armata-SBL558 Eugeissona-tristis-MSL52-S20 Calamus-ornatus-BKL092 Plectocomiopsis-mira-RBL352 Mauritia-flexuosa-MSL70 Calamus-ursinus-BKL042 Calamus-essigii-RBL130 Korthalsia-rostrata-RBL177 Laccosperma-opacum-MSL50 Pigafetta-filaris-MSL77-S28 Calamus-erectus-BKL037 Calamus-ciliaris-BKL028 Calamus-symphysipus-RBL089 Calamus-thysanolepis-BKL160 Eremospatha-laurentii-BKL002 Calamus-pedicellatus-RBL031 Eremospatha-wendlandiana-BKL087 Calamus-conirostris-RBL214 Calamus-zollingeri-RBL145 Calamus-aruensis-RBL231 Salacca-zalacca-BKL003 Calamus-subinermis-BKL018 Calamus-wailong-RBL021 Calamus-arborescens-RBL069 Calamus-harmandii-BKL025 Oncocalamus-tuleyi-BKL016 Raphia-sudanica-BKL044 Calamus-vitiensis-RBL211 Calamus-blumei-RBL084 Calamus-pseudoconcolor-RBL081 Calamus-usitatus-BKL019 Calamus-koordersianus-BKL031 Calamus-acanthophyllus-BKL024 Calamus-dumetosus-RBL215 Calamus-melanochaetes-RBL238 Calamus-pogonacanthus-BKL048 Calamus-rhabdocladus-RBL032 Korthalsia-jala-BKL014 Eugeissona-utilis-BKL081 Calamus-godefroyi-BKL128 Salacca-affinis-BKL055 Calamus-oxleyanus-BKL033 Calamus-compsostachys-RBL028 Laccosperma-secundiflorum-RBL066 Calamus-castaneus-RBL220 Korthalsia-rigida-RBL038 Calamus-acanthochlamys-RBL017 Pigafetta-elata-BKL084 Calamus-sabut-RBL024 Salacca-secunda-BKL077 Metroxylon-salomonense-BKL060 Calamus-deerratus-RBL242 Korthalsia-robusta-BKL083 -f fasta -d dna -i alignment_genes_concatenated.out```
