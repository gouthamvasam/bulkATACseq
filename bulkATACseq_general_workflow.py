# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 12:09:46 2024

@author: gouthamvasam
"""

# Import necessary libraries
import subprocess

# Quality Control of Raw Reads using FastQC
def run_fastqc(fastq_files, output_dir):
    fastqc_cmd = f"fastqc {' '.join(fastq_files)} --outdir={output_dir}"
    subprocess.run(fastqc_cmd, shell=True, check=True)

# Read Alignment using Bowtie2
def align_reads(fastq_files, output_bam, reference_genome):
    # Align reads
    bowtie2_cmd = f"bowtie2 -x {reference_genome} -U {' '.join(fastq_files)} | samtools view -bS - > {output_bam}"
    subprocess.run(bowtie2_cmd, shell=True, check=True)
    # Sort BAM file
    sorted_bam = output_bam.replace('.bam', '.sorted.bam')
    sort_cmd = f"samtools sort {output_bam} -o {sorted_bam}"
    subprocess.run(sort_cmd, shell=True, check=True)
    # Index BAM file
    index_cmd = f"samtools index {sorted_bam}"
    subprocess.run(index_cmd, shell=True, check=True)
    return sorted_bam

# Removal of PCR Duplicates using Picard
def remove_duplicates(input_bam, output_bam):
    picard_cmd = f"java -jar picard.jar MarkDuplicates I={input_bam} O={output_bam} M=marked_dup_metrics.txt REMOVE_DUPLICATES=true"
    subprocess.run(picard_cmd, shell=True, check=True)
    
def filter_mitochondrial(input_bam, output_bam):
    # Create a file with the names of all chromosomes except the mitochondrial chromosome
    exclude_mito_cmd = f"samtools idxstats {input_bam} | grep -v '^chrM' | cut -f1 > exclude_chrM.txt"
    subprocess.run(exclude_mito_cmd, shell=True, check=True)
    # Filter out mitochondrial reads using the file
    filter_cmd = f"samtools view -b {input_bam} -o temp_no_mito.bam -L exclude_chrM.txt"
    subprocess.run(filter_cmd, shell=True, check=True)
    # Sort and index the BAM file after mitochondrial reads are filtered out
    sort_cmd = f"samtools sort temp_no_mito.bam -o {output_bam}"
    subprocess.run(sort_cmd, shell=True, check=True)
    index_cmd = f"samtools index {output_bam}"
    subprocess.run(index_cmd, shell=True, check=True)
    # Clean up temporary files
    subprocess.run("rm temp_no_mito.bam exclude_chrM.txt", shell=True, check=True)

# Peak Calling using MACS2
def call_peaks(input_bam, output_dir):
    macs2_cmd = f"macs2 callpeak -t {input_bam} --outdir {output_dir} --name peaks --format BAMPE --nomodel --nolambda --keep-dup all --call-summits"
    subprocess.run(macs2_cmd, shell=True, check=True)

# Create a tag directory for downstream analysis using HOMER
def create_tag_directory(input_bam, tag_directory):
    make_tag_directory_cmd = f"makeTagDirectory {tag_directory} {input_bam}"
    subprocess.run(make_tag_directory_cmd, shell=True, check=True)

# Perform motif analysis using HOMER
def find_motifs_genome(tag_directory, output_dir, genome):
    find_motifs_cmd = f"findMotifsGenome.pl {tag_directory} {genome} {output_dir} -size 200 -mask"
    subprocess.run(find_motifs_cmd, shell=True, check=True)

# Annotate peaks using HOMER
def annotate_peaks(peaks_file, genome, output_file):
    annotate_peaks_cmd = f"annotatePeaks.pl {peaks_file} {genome} > {output_file}"
    subprocess.run(annotate_peaks_cmd, shell=True, check=True)

# Create bigWig files for visualization using deepTools
def bam_to_bigwig(input_bam, genome_size, output_bw):
    bam_coverage_cmd = f"bamCoverage -b {input_bam} -o {output_bw} --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize {genome_size}"
    subprocess.run(bam_coverage_cmd, shell=True, check=True)

# Create a coverage matrix for the heatmap
def create_coverage_matrix(bigwig_files, regions_file, output_matrix):
    compute_matrix_cmd = f"computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -R {regions_file} -S {' '.join(bigwig_files)} -out {output_matrix}"
    subprocess.run(compute_matrix_cmd, shell=True, check=True)

# Plot a heatmap using the coverage matrix
def plot_heatmap(matrix_file, output_heatmap):
    plot_heatmap_cmd = f"plotHeatmap -m {matrix_file} -out {output_heatmap}"
    subprocess.run(plot_heatmap_cmd, shell=True, check=True)

# Plot a profile plot using the coverage matrix
def plot_profile(matrix_file, output_profile):
    plot_profile_cmd = f"plotProfile -m {matrix_file} -out {output_profile}"
    subprocess.run(plot_profile_cmd, shell=True, check=True)
    
# Main function to run the workflow
def main():
    # Example usage:
    fastq_files = ['sample1.fastq', 'sample2.fastq']
    output_dir = 'fastqc_results'
    reference_genome = 'genome_index'
    genome_size = 'hs'  # Replace with the actual effective genome size for your organism
    sorted_bam = 'sorted_reads.bam'
    dedup_bam = 'dedup_reads.bam'
    no_mito_bam = 'no_mito_reads.bam'
    macs2_output_dir = 'macs2_output'
    tag_directory = 'ATAC_tagDir'
    motif_output_dir = 'HOMER_motifs'
    peaks_file = f"{macs2_output_dir}/peaks.narrowPeak"
    annotated_peaks_file = 'annotated_peaks.txt'
    bigwig_file = 'reads.bw'
    regions_file = peaks_file  # Using the output from MACS2 as the regions file
    coverage_matrix = 'coverage_matrix.gz'
    heatmap_output = 'heatmap.png'
    profile_output = 'profile.png'

    # Run the workflow
    run_fastqc(fastq_files, output_dir)
    sorted_bam = align_reads(fastq_files, sorted_bam, reference_genome)
    remove_duplicates(sorted_bam, dedup_bam)
    filter_mitochondrial(dedup_bam, no_mito_bam)
    call_peaks(no_mito_bam, macs2_output_dir)
    create_tag_directory(no_mito_bam, tag_directory)
    find_motifs_genome(tag_directory, motif_output_dir, reference_genome)
    annotate_peaks(peaks_file, reference_genome, annotated_peaks_file)
    bam_to_bigwig(no_mito_bam, genome_size, bigwig_file)
    create_coverage_matrix([bigwig_file], regions_file, coverage_matrix)
    plot_heatmap(coverage_matrix, heatmap_output)
    plot_profile(coverage_matrix, profile_output)

if __name__ == "__main__":
    main()