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
