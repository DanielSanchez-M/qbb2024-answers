## Q1
- grep "protein_coding" hg38-gene-metadata-feature.tsv | cut -f 7 | uniq -c             
- 19618 protein_coding

or

- cut -f 7 hg38-gene-metadata-feature.tsv | sort | uniq -c
# 19618 Protein Coding Genes
# Other Bio_type = miRNA = 1833
# a single small miRNA gene has the potential to regulate 100s of other genes during transcription. Can be very critical to understanding development and transcription processes.