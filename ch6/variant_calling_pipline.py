"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
###########################################
This Python script implements a complete germline variant calling pipeline using the GATK Best Practices 
workflow. It performs the following steps: reference genome download and indexing, quality control, 
read trimming, alignment, duplicate marking, base quality score recalibration, variant calling using 
HaplotypeCaller, joint genotyping, and variant filtration.

The script leverages multiple tools including BWA, SAMtools, Picard, GATK, Trimmomatic, and FastQC.
It is designed to automate the analysis of multiple samples using metadata stored in a CSV file, 
with results organized per sample and final cohort VCFs generated at the end.

Required Packages:
- Python: pandas, subprocess, os
- External tools (must be installed and accessible in PATH):
  - FastQC
  - Trimmomatic
  - BWA
  - SAMtools
  - GATK
  - Picard
  - wget

Input Files:
- Reference genome URL (e.g., hg38 FASTA from UCSC): downloaded and indexed automatically
- Known variants VCF file (e.g., hapmap_3.3.hg38.vcf.gz): required for base recalibration
- Metadata CSV file: meta/metadata.csv with a 'runID' column for sample identifiers
- Paired-end FASTQ files: data/raw/sampleID_1.fastq.gz and sampleID_2.fastq.gz for each sample

Output Files:
- Reference directory: indexed genome FASTA, dict, and BWA index (reference/)
- Quality control reports: results/sampleID/qc/
- Trimmed reads: results/sampleID/trimmed/
- Aligned and processed BAMs: results/sampleID/aligned/ and results/sampleID/processed/
- GVCFs for each sample: results/sampleID/gvcf/
- Joint genotyped VCF: results/cohort/cohort.vcf.gz
- Filtered final VCF: results/cohort/filtered_variants.vcf.gz

Note:
Make sure that all external tools are installed, and `known_vcfs/hapmap_3.3.hg38.vcf.gz` exists before running.
"""

import subprocess
import os
import pandas as pd

def run_command(command):
    print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True)

def download_and_index_reference(reference_url, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    compressed_path = os.path.join(output_dir, "reference.fasta.gz")
    ref_fasta = os.path.join(output_dir, "reference.fasta")
    ref_dict = ref_fasta.replace(".fasta", ".dict")
    run_command(f"wget -O {compressed_path} {reference_url}")
    run_command(f"gunzip -c {compressed_path} > {ref_fasta}")
    run_command(f"samtools faidx {ref_fasta}")
    # Delete .dict if it exists
    if os.path.exists(ref_dict):
        os.remove(ref_dict)
    run_command(
     f"gatk CreateSequenceDictionary "
     f"-R {ref_fasta} "
     f"-O {ref_dict}"
    )
    run_command(f"bwa index {ref_fasta}")
    return ref_fasta

def quality_control(fastq_files, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for fq in fastq_files:
        run_command(f"fastqc {fq} -o {output_dir}")

def trim_reads(sample_name, r1, r2, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    out_r1 = os.path.join(output_dir, f"{sample_name}_R1_trimmed.fastq.gz")
    out_r2 = os.path.join(output_dir, f"{sample_name}_R2_trimmed.fastq.gz")
    run_command(
        f"trimmomatic PE {r1} {r2} {out_r1} /dev/null {out_r2} /dev/null "
        f"SLIDINGWINDOW:4:20 MINLEN:50"
    )
    return out_r1, out_r2

def align_reads(sample_name, ref_genome, r1_trimmed, r2_trimmed, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    sam_file = os.path.join(output_dir, f"{sample_name}.sam")
    bam_file = os.path.join(output_dir, f"{sample_name}.bam")
    rg_bam_file = os.path.join(output_dir, f"{sample_name}_rg.bam")
    # Step 1: Align without read group
    run_command(
     f"bwa mem {ref_genome} "
     f"{r1_trimmed} {r2_trimmed} "
     f"> {sam_file}"
    )
    # Step 2: Convert to BAM
    run_command(f"samtools view -Sb {sam_file} > {bam_file}")
    # Step 3: Add proper read group
    run_command(
    (
        f"picard AddOrReplaceReadGroups "
        f"I={bam_file} "
        f"O={rg_bam_file} "
        f"RGID={sample_name} "
        f"RGLB=lib1 "
        f"RGPL=ILLUMINA "
        f"RGPU=unit1 "
        f"RGSM={sample_name} "
        f"VALIDATION_STRINGENCY=LENIENT"
       )
    )  
    return rg_bam_file

def convert_sort_mark_duplicates(sample_name, bam_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    sorted_bam = os.path.join(output_dir, f"{sample_name}_sorted.bam")
    dedup_bam = os.path.join(output_dir, f"{sample_name}_dedup.bam")
    metrics_file = os.path.join(output_dir, f"{sample_name}_metrics.txt")
    run_command(f"samtools sort {bam_file} -o {sorted_bam}")
    run_command(
        f"gatk MarkDuplicates "
        f"-I {sorted_bam} "
        f"-O {dedup_bam} "
        f"-M {metrics_file}"
    )
    return dedup_bam

def base_recalibration(sample_name, dedup_bam, ref_genome, known_sites_vcf, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    recal_table = os.path.join(output_dir, f"{sample_name}_recal_data.table")
    recal_bam = os.path.join(output_dir, f"{sample_name}_recal.bam")
    run_command(
        f"gatk BaseRecalibrator "
        f"-I {dedup_bam} "
        f"-R {ref_genome} "
        f"--known-sites {known_sites_vcf} "
        f"-O {recal_table}"
    )
    run_command(
        f"gatk ApplyBQSR "
        f"-I {dedup_bam} "
        f"-R {ref_genome} "
        f"--bqsr-recal-file {recal_table} "
        f"-O {recal_bam}"
    )
    return recal_bam

def haplotype_caller(sample_name, recal_bam, ref_genome, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    gvcf = os.path.join(output_dir, f"{sample_name}.g.vcf.gz")
    run_command(
        f"gatk HaplotypeCaller "
        f"-R {ref_genome} "
        f"-I {recal_bam} "
        f"-O {gvcf} "
        f"-ERC GVCF"
    )
    return gvcf

def joint_genotyping(gvcf_list, ref_genome, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    merged_gvcf = os.path.join(output_dir, "cohort.g.vcf.gz")
    vcf_output = os.path.join(output_dir, "cohort.vcf.gz")
    inputs = ' '.join([f"--variant {gvcf}" for gvcf in gvcf_list])
    run_command(
        f"gatk CombineGVCFs -R {ref_genome} {inputs} -O {merged_gvcf}"
    )
    run_command(
        f"gatk GenotypeGVCFs -R {ref_genome} -V {merged_gvcf} -O {vcf_output}"
    )
    return vcf_output

def filter_variants(vcf_file, ref_genome, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    filtered_vcf = os.path.join(output_dir, "filtered_variants.vcf.gz")
    run_command(
        f"gatk VariantFiltration "
        f"-R {ref_genome} "
        f"-V {vcf_file} "
        f"-O {filtered_vcf} "
        f"--filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0' "
        f"--filter-name 'GATK_Hard_Filter'"
    )
    return filtered_vcf

# Run the pipeline
if __name__ == '__main__':
    reference_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    reference_dir = "reference"
    ref_genome = download_and_index_reference(reference_url, reference_dir)
    known_sites_vcf = "known_vcfs/hapmap_3.3.hg38.vcf.gz"

    metadata = pd.read_csv("meta/metadata.csv")

    all_gvcfs = []
    for _, row in metadata.iterrows():
        sample = row['runID']
        r1 = f"data/raw/{sample}_1.fastq.gz"
        r2 = f"data/raw/{sample}_2.fastq.gz"
        qc_dir = f"results/{sample}/qc"
        trimmed_dir = f"results/{sample}/trimmed"
        aligned_dir = f"results/{sample}/aligned"
        processed_dir = f"results/{sample}/processed"
        gvcf_dir = f"results/{sample}/gvcf"
        quality_control([r1, r2], qc_dir)
        r1_trimmed, r2_trimmed = trim_reads(sample, r1, r2, trimmed_dir)
        bam_with_rg = align_reads(sample, ref_genome, r1_trimmed, r2_trimmed, aligned_dir)
        dedup_bam = convert_sort_mark_duplicates(sample, bam_with_rg, processed_dir)
        recal_bam = base_recalibration(sample, dedup_bam, ref_genome, known_sites_vcf, processed_dir)
        gvcf = haplotype_caller(sample, recal_bam, ref_genome, gvcf_dir)
        all_gvcfs.append(gvcf)
    joint_output_dir = "results/cohort"
    final_vcf = joint_genotyping(all_gvcfs, ref_genome, joint_output_dir)
    filtered_vcf = filter_variants(final_vcf, ref_genome, joint_output_dir)
    print(f"Pipeline complete. Final annotated VCF: {annotated_vcf}")

