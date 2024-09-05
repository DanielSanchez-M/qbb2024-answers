# Q1
- (qb24) cmdb@QuantBio-10 d2_exercises %grep "protein_coding" hg38-gene-metadata-feature.tsv | cut -f 7 | uniq -c             
    -  19618 protein_coding
- (qb24) cmdb@QuantBio-10 d2_exercises % cut -f 7 hg38-gene-metadata-feature.tsv | sort | uniq -c

    - 19618 Protein Coding Genes
    - Other Bio_type = miRNA = 1833
    - a single small miRNA gene has the potential to regulate 100s of other genes during transcription. Can be very critical to understanding development and transcription processes.

# Q2
- (qb24) cmdb@QuantBio-10 d2_exercises % cut -f 1 hg38-gene-metadata-go.tsv | uniq -c | sort -n
    - Most GO terms =  273 ENSG00000168036

- (qb24) cmdb@QuantBio-10 d2_exercises % grep "ENSG00000168036" hg38-gene-metadata-go.tsv | sort -k 3 > gene_ENSG00000168036.tsv
    - This gene is a realy important developmental gene given by the number of processes they are involved in, some of which are involved embryonic morphogenesis. Ex: Protein binding, embryonic brain development, embryonic foregut morphogenesis, etc

# Q1
- (qb24) cmdb@QuantBio-10 d2_exercises % grep -e "IG_._gene" gene.gtf | cut -f 1 | uniq -c | sort 
   - 1 chr21
   6 chr16
  16 chr15
  48 chr22
  52 chr2
  91 chr14

- (qb24) cmdb@QuantBio-10 d2_exercises % grep -e "IG_._pseudogene" -e "IG_pseudogene" gene.gtf | cut -f 1 | uniq -c | sort
    - 1 chr1
   1 chr10
   1 chr18
   1 chr8
   5 chr9
   6 chr15
   8 chr16
  45 chr2
  48 chr22
  84 chr14
    - The distribution of IG pseudogenes are spread across more chromosomes than IG genes. While chromosomes 22, 14, 15, 16, 2 all have both IG genes and IG pseudogenes present.

# Q2
- grep pseudogene gene.gtf will not be specific to the IG genes and instead gather the data from all those classified as pseudogene. Including the TR genes, non-coding genes,etc.
- A better way to search for the specific pseudogene under gene_type would be to use the grep function and specify specifically which terms should be included and which ones should be excluded. The terms overlap would be added to be excluded.
- (qb24) cmdb@QuantBio-10 d2_exercises % grep -F "pseudogene" genes.gtf | grep  -v "overlap" | less -S

# Q3
- (base) cmdb@QuantBio-10 d2_exercises % cut -f1,4,5,14 gene-tabs.gtf > genes.bed
