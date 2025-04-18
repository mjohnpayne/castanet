# Adjust channel priority on conda env
conda config --add channels bioconda
conda config --add channels conda-forge

# Install kraken2 via conda, for removing human reads
conda install -y -c bioconda kraken2=2.1.3

# Install samtools
conda install -y "samtools>=1.10"

# Install mafft
conda install -y -c bioconda mafft=7.520

# # Download pre-built kraken2 database with human genome only
mkdir kraken2_human_db
curl -L -o kraken2_human_db/kraken2_human_db.tar.gz https://ndownloader.figshare.com/files/23567780
tar -xzvf kraken2_human_db/kraken2_human_db.tar.gz

# Download trimmomatic and extract
conda install -y trimmomatic=0.39

# Install mappers
conda install -y bwa-mem2=2.2.1
conda install -y bowtie2=2.5.4

# Install viral consensus tool
conda install -y viral_consensus=0.0.5
