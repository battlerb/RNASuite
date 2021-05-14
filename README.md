# rna-seq-pipeline

#### Introduction
pipeline for analysis RNA-seq data.

Version: 0.2.0(BETA)

#### Pipeline
This pipeline include multi submodules: quality contral, alignment, call count, DEG analysis and enrichment.

#### Prepare
    git clone https://github.com/battlerb/RNASuite.git && cd RNASuite && Rscript rna-seq-pipeline/install_packages.R

#### Guide
Please download reference genome and annotation gtf file first and move them to `/db/Genome/GRC*38/` folder, recommand Ensembl release. 

Human： [fasta](ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) [gtf](ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz)

Mouse： [fasta](ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz) [gtf](ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.chr.gtf.gz)

Pig：[fasta](ftp://ftp.ensembl.org/pub/release-101/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz) [gtf](ftp://ftp.ensembl.org/pub/release-101/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.101.chr.gtf.gz)

#### Run

    python3 wrapper.py -l sampleList.txt -p 8 -o output -s hsa -i fastq_folder/ -c compare.csv