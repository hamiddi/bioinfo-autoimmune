#Author: Hamid D. Ismail, Ph.D.
#Book: Bioinformatics of Autoimmune Diseases

#This shell script renames simplified paired-end FASTQ files (e.g., sample_1.fastq.gz and sample_2.fastq.gz) 
#into the full Illumina CASAVA 1.8+ format (e.g., sample_S1_L001_R1_001.fastq.gz and sample_S1_L001_R2_001.fastq.gz). 
#This renaming is often required by pipelines or tools expecting CASAVA-style filenames, such as Illumina's bcl2fastq 

cd data/raw
for file in *_1.fastq.gz; do
    base=$(echo $file | sed 's/_1\.fastq\.gz//')
    mv ${base}_1.fastq.gz ${base}_S1_L001_R1_001.fastq.gz
    mv ${base}_2.fastq.gz ${base}_S1_L001_R2_001.fastq.gz
done
