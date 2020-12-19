# PhyloFrame
Included here are scripts to replicate phylogenomic analyses in:

Kuhnhäuser et al. (2021), A robust phylogenomic framework for the calamoid palms, Molecular Phylogenetics and Evolution.

Raw sequence data are deposited in the European Nucleotide Archive of the European Bioinformatics Institute (https://www.ebi.ac.uk/ena) under project number PRJEB40689. The PhyloPalm target file and alignments, gene trees, and species trees resulting from the analyses detailed here are deposited on Zenodo (https://doi.org/10.5281/zenodo.4359280).

Notes:
- Exemplary workflow for taxon "Calamoid1" and exon "Gene1"
- Identical parameters used for exons and supercontigs

## 1. Quality control of sequence data
```fastqc Calamoid1_R*_001.fastq --extract```
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

## 8. Gene trees
### Without model testing
```raxmlHPC-PTHREADS -T 3 -m GTRGAMMA -f a -p 12345 -x 12345 -# 1000 -k -s Gene1_aligned_trimmed.fasta -n Gene1.tree```
### With model testing
```iqtree -s Gene1_aligned_trimmed.fasta -m MFP -T 4 -B 1000```
- ```-m MFP``` performs model testing followed by tree search

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
#### Without partitioning
1) Concatenate individual gene alignments\
```AMAS concat -i Gene*_aligned_trimmed.fasta -f fasta -d dna -c 1```
2) Build concatenated species tree\
```raxmlHPC-PTHREADS -T 3 -m GTRGAMMA -f a -p 12345 -x 12345 -# autoMRE -k -s concatenated.out -n raxml_concatenated.tree```
- ```autoMRE``` argument for extended majority-rule criterion for stopping bootstrapping after a sufficient number of replicates had been sampled 
3) Root the tree\
```pxrr -t raxml_concatenated.tree -g Nypa-fructicans-MSL30-S32,Kerriodoxa-elegans-MSL76,Asterogyne-martiana-SBL226,Ceroxylon-quindiuense-MSL17 -s > raxml_concatenated_rooted.tree```

#### With partitioning
1) Concatenate individual gene alignments\
```AMAS concat -i Gene*_aligned_trimmed.fasta -f fasta -d dna -c 1```
2) Perform PartitionFinder analysis\
```partitionfinder-2.1.1/PartitionFinder.py partitioned_analysis --raxml```\
Directory ```partitioned_analysis``` contains:\
- concatenated alignment in phylyp format
- configuration file ```partition_finder.cfg``` with partitions file from ```AMAS concat``` output in DATA BLOCKS:
```
## ALIGNMENT FILE ##
alignment = concatenated.phy;

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##
models = GTR+G;

# MODEL SELECCTION: AIC | AICc | BIC #
model_selection = aicc;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]
p1_EGU105032175_supercontig_aligned_trimmed_sites_seq = 1-675;
...
p945_HEY989_supercontig_aligned_trimmed_sites_seq = 1942327-1943739;

## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##
[schemes]
search = rcluster;
```
- Copy RAxML partitions from output file "best_scheme.txt" to file "partitions_partitionfinder.txt" as input for RAxML analyses
3) Partitioned RAxML analysis 
```raxmlHPC-PTHREADS -T 22 -m GTRGAMMA -f a -p 12345 -x 12345 -# autoMRE -k -s concatenated.out -n raxml_supercontigs_concatenated_partitioned_partitionfinder_autoMRE.tree -q partitions_partitionfinder.txt```

## 10. Comparative analyses
### Tree comparisons
#### Load R libraries:\
```library(ape)```\
```library(phangorn)```\
#### Read in trees:\
```coal_e <- read.tree("Calamoideae_coalescence_exons_rooted.tree")```\
```coal_em <- read.tree("Calamoideae_coalescence_exons_models_rooted.tree")```\
```coal_s <- read.tree("Calamoideae_coalescence_supercontigs_rooted.tree")```\
```coal_sm <- read.tree("Calamoideae_coalescence_supercontigs_models_rooted.tree")```\
```concat_e <- read.tree("Calamoideae_concatenation_exons_rooted.tree")```\
```concat_ep <- read.tree("Calamoideae_concatenation_exons_partitioned_rooted.tree")```\
```concat_s <- read.tree("Calamoideae_concatenation_supercontigs_rooted.tree")```\
```concat_sp <- read.tree("Calamoideae_concatenation_supercontigs_partitioned_rooted.tree")```

#### Strict consensus of species trees
```cons <- consensus(coal_e, coal_em, coal_s, coal_sm, concat_e, concat_ep, concat_s, concat_sp)```\
```plot.phylo(cons, cex = 0.7)```\
```write.tree(cons, file = "consensus_all.tree")```

#### Normalized Robinson Foulds distances
```all_trees <- c(coal_e, coal_em, coal_s, coal_sm, concat_e, concat_ep, concat_s, concat_sp)```\
```RF.dist(all_trees, check.labels = TRUE, normalize = TRUE, rooted = TRUE)```

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
- For outgroup only: first remove ingroup (by listing all ingroup taxa; analogous to "ingroup only" analysis)
