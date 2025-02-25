"""
Author: YF Sun
Date: 2024-06-05
Description: CRISPR Sequence Mapping and Deletion Analysis Pipeline
"""

import argparse
import os
import subprocess
from pathlib import Path  # Better path handling
from Bio.Seq import Seq

# Constants (should be moved to config file or command-line arguments)
TOOL_PATHS = {
    'cutadapt': "/path/to/cutadapt",
    'seqkit': "/path/to/seqkit",
    'star': "/path/to/STAR",
    'python': "/path/to/python"
}

TOOL_PATHS = {
    'cutadapt': "/share/public4/data/sunyf/miniconda3/envs/python38/bin/cutadapt",
    'seqkit': "/share/public4/data/sunyf/miniconda3/envs/python38/bin/seqkit", 
    'star': "/share/public4/data/sunyf/miniconda3/envs/python38/bin/STAR",
    'python': "/share/public4/data/sunyf/miniconda3/envs/python38/bin/python3.8"
}

def create_sample_directories(work_dir: str, sample: str) -> Path:
    """Create directory structure for sample processing.
    
    Args:
        work_dir: Parent working directory
        sample: Sample identifier
    
    Returns:
        Path object to sample directory
    """
    sample_dir = Path(work_dir) / sample
    temp_dir = sample_dir / "temp"
    
    sample_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(exist_ok=True)
    
    return sample_dir

def process_adapters(input_reads: str, output_dir: Path, 
                    remove_seq: str = None, retain_seq: str = None) -> str:
    """Process sequencing reads with adapter removal and filtering.
    
    Args:
        input_reads: Path to input FASTQ file
        output_dir: Output directory path
        remove_seq: Sequence pattern to remove
        retain_seq: Sequence pattern to retain
    
    Returns:
        Path to processed FASTQ file
    """
    processed_reads = input_reads
    
    # Sequence filtering pipeline
    if remove_seq:
        output_rm = output_dir / "reads.seqkit_rm.fastq"
        cmd = f"{TOOL_PATHS['seqkit']} grep -s -i -j 10 -P -m 1 -p {remove_seq} -v -o {output_rm} {processed_reads}"
        subprocess.run(cmd, shell=True, check=True)
        processed_reads = output_rm

    if retain_seq:
        output_rt = output_dir / "reads.seqkit_rt.fastq"
        cmd = f"{TOOL_PATHS['seqkit']} grep -s -i -j 10 -P -m 1 -p {retain_seq} -o {output_rt} {processed_reads}"
        subprocess.run(cmd, shell=True, check=True)
        processed_reads = output_rt

    # Adapter trimming
    output_cut = output_dir / "read.cutadapt.fastq"
    cmd = (
        f"{TOOL_PATHS['cutadapt']} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA "
        f"-j 10 --quiet --max-n 0 -e 0.1 -q 20 -m 20 --trim-n -o {output_cut} {processed_reads}"
    )
    subprocess.run(cmd, shell=True, check=True)
    
    return output_cut

def build_star_index(reference_seq: str, output_dir: Path) -> None:
    """Build STAR genome index from reference sequence.
    
    Args:
        reference_seq: Reference DNA sequence
        output_dir: Directory for index files
    """
    index_file = output_dir / "index.fa"
    with open(index_file, "w") as f:
        f.write(f">raw\n{reference_seq}\n")
    
    cmd = (
        f"{TOOL_PATHS['star']} --runMode genomeGenerate "
        f"--runThreadN 10 --genomeDir {output_dir} "
        f"--genomeSAindexNbases 3 --genomeFastaFiles {index_file} "
        f"1> {output_dir}/std_out_indexing"
    )
    subprocess.run(cmd, shell=True, check=True)

def run_star_alignment(processed_reads: str, index_dir: Path, 
                      output_dir: Path, num_threads: int = 5) -> None:
    """Run STAR alignment with optimized parameters.
    
    Args:
        processed_reads: Path to processed FASTQ file
        index_dir: Path to STAR index directory
        output_dir: Output directory for alignment results
        num_threads: Number of threads to use
    """
    cmd = (
        f"{TOOL_PATHS['star']} --genomeDir {index_dir} "
        f"--runThreadN {num_threads} --readFilesIn {processed_reads} "
        f"--outFileNamePrefix {output_dir}/EndToEnd_to_dc_ "
        f"--outReadsUnmapped Fastx --outSAMprimaryFlag AllBestScore "
        f"--alignEndsType EndToEnd --scoreDelOpen 0 --scoreDelBase 0 "
        f"--clip5pNbases 10 --alignIntronMin 100 1> {output_dir}/std_out_mapping"
    )
    subprocess.run(cmd, shell=True, check=True)

def analyze_deletions(bam_file: Path, reference_seq: str, guide_rna: str,
                     sample_name: str, output_dir: Path) -> None:
    """Analyze deletion patterns from aligned reads.
    
    Args:
        bam_file: Path to sorted BAM file
        reference_seq: Reference DNA sequence
        guide_rna: Guide RNA sequence
        sample_name: Sample identifier
        output_dir: Output directory for results
    """
    cmd = (
        f"{TOOL_PATHS['python']} /path/to/Call_deletion_version4.2.py "
        f"--input_bam {bam_file} --gRNA_seq {guide_rna} "
        f"--raw_seq {reference_seq} --sample_name {sample_name} "
        f"--output_csv {output_dir}/df_level_version4.2.csv"
    )
    subprocess.run(cmd, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description='CRISPR Sequence Analysis Pipeline')
    
    # Required arguments
    parser.add_argument('--sample_name', required=True, help='Sample identifier')
    parser.add_argument('--read1', required=True, help='Path to read 1 FASTQ')
    parser.add_argument('--read2', required=True, help='Path to read 2 FASTQ')
    parser.add_argument('--raw_seq', required=True, help='Reference DNA sequence')
    parser.add_argument('--gRNA', required=True, help='Guide RNA sequence')
    
    # Optional arguments
    parser.add_argument('--work_dir', default='.', help='Working directory')
    parser.add_argument('--threads', type=int, default=5, help='Number of threads')
    parser.add_argument('--steps', type=int, choices=[1,12,123,3,4], 
                       help='Processing steps to execute')
    
    args = parser.parse_args()
    
    # Setup directory structure
    sample_dir = create_sample_directories(args.work_dir, args.sample_name)
    
    # Execute processing steps
    if args.steps in [12, 123]:
        processed_reads = process_adapters(
            args.read1, sample_dir, 
            args.sequence_remove, args.sequence_retain
        )
    
    if args.steps in [1, 123]:
        build_star_index(args.raw_seq.upper(), sample_dir)
    
    if args.steps == 123:
        run_star_alignment(processed_reads, sample_dir, sample_dir, args.threads)
        bam_file = sample_dir / "EndToEnd_to_dc_Aligned.out.sorted.bam"
        analyze_deletions(bam_file, args.raw_seq.upper(), 
                         args.gRNA.upper(), args.sample_name, sample_dir)

if __name__ == "__main__":
    main()