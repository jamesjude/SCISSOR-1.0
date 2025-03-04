"""
Author: YF Sun
Date: 2024-06-05
Description: CRISPR Sequence Mapping and Deletion Analysis Pipeline

"""

import os
import re
import csv
import pysam
import subprocess
import sys
from pathlib import Path
from typing import Optional, Dict, List, Tuple
from Bio.Seq import Seq
from collections import defaultdict, Counter
from tqdm import tqdm
import argparse

# Configuration for external tool paths
TOOL_PATHS = {
    'cutadapt': "/share/public4/data/sunyf/miniconda3/envs/python38/bin/cutadapt",
    'seqkit': "/share/public4/data/sunyf/miniconda3/envs/python38/bin/seqkit",
    'star': "/share/public4/data/sunyf/miniconda3/envs/python38/bin/STAR",
    'python': "/share/public4/data/sunyf/miniconda3/envs/python38/bin/python3.8"
}

class PipelineError(Exception):
    """Custom exception for pipeline failures"""
    pass

def validate_tools(tool_paths: Dict[str, str]) -> None:
    """Validate existence of required external tools
    
    Args:
        tool_paths: Dictionary mapping tool names to their paths
        
    Raises:
        PipelineError: If any required tool is not found
    """
    for tool_name, path in tool_paths.items():
        if not Path(path).exists():
            raise PipelineError(f"Required tool not found: {tool_name} at {path}")


def create_sample_directories(work_dir: Path, sample: str) -> Path:
    """Create directory structure for sample processing
    
    Args:
        work_dir: Parent working directory path
        sample: Unique sample identifier
        
    Returns:
        Path: Created sample directory path
        
    Raises:
        PipelineError: If directory creation fails
    """
    try:
        sample_dir = work_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        (sample_dir / "temp").mkdir(exist_ok=True)
        return sample_dir
    except OSError as e:
        raise PipelineError(f"Directory creation failed: {str(e)}")

def process_reads(
    input_reads: Path,
    output_dir: Path,
    remove_seq: Optional[str] = None,
    retain_seq: Optional[str] = None
) -> Path:
    """Process sequencing reads with adapter removal and filtering
    
    Args:
        input_reads: Path to input FASTQ file
        output_dir: Output directory path
        remove_seq: Sequence pattern to remove (optional)
        retain_seq: Sequence pattern to retain (optional)
        
    Returns:
        Path: Path to processed FASTQ file
        
    Raises:
        PipelineError: If processing fails at any stage
    """
    current_reads = input_reads
    
    try:
        # Sequence filtering operations
        operations = []
        if remove_seq:
            operations.append(('remove', remove_seq, '-v'))
        if retain_seq:
            operations.append(('retain', retain_seq, ''))
            
        for op_name, seq, flag in operations:
            output_file = output_dir / f"reads.seqkit_{op_name}.fastq"
            cmd = [
                TOOL_PATHS['seqkit'], 'grep',
                '-s', '-i', '-j', '10', '-P', '-m', '1',
                '-p', seq, flag, '-o', str(output_file),
                str(current_reads)
            ]
            subprocess.run(cmd, check=True)
            current_reads = output_file

        # Adapter trimming with cutadapt
        output_cut = output_dir / "read.cutadapt.fastq"
        cmd = [
            TOOL_PATHS['cutadapt'],
            '-a', 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
            '-j', '10', '--quiet', '--max-n', '0',
            '-e', '0.1', '-q', '20', '-m', '20',
            '--trim-n', '-o', str(output_cut),
            str(current_reads)
        ]
        subprocess.run(cmd, check=True)
        
        return output_cut
        
    except subprocess.CalledProcessError as e:
        raise PipelineError(f"Read processing failed at step: {e.cmd}")

def build_star_index(reference_seq: str, index_dir: Path) -> None:
    """Build STAR genome index from reference sequence
    
    Args:
        reference_seq: Reference DNA sequence
        index_dir: Directory for index files
        
    Raises:
        PipelineError: If index creation fails
    """
    try:
        index_file = index_dir / "index.fa"
        with open(index_file, "w") as f:
            f.write(f">raw\n{reference_seq.upper()}\n")
            
        cmd = [
            TOOL_PATHS['star'], '--runMode', 'genomeGenerate',
            '--runThreadN', '10', '--genomeDir', str(index_dir),
            '--genomeSAindexNbases', '3', '--genomeFastaFiles', str(index_file)
        ]
        subprocess.run(cmd, check=True, capture_output=True)
        
    except (subprocess.CalledProcessError, OSError) as e:
        raise PipelineError(f"Index creation failed: {str(e)}")

def run_star_alignment(
    processed_reads: Path,
    index_dir: Path,
    output_dir: Path,
    num_threads: int = 5
) -> Path:
    """Run STAR alignment with optimized parameters
    
    Args:
        processed_reads: Path to processed FASTQ file
        index_dir: Path to STAR index directory
        output_dir: Output directory for alignment results
        num_threads: Number of threads to use
        
    Returns:
        Path: Path to generated BAM file
        
    Raises:
        PipelineError: If alignment fails
    """
    try:
        cmd = [
            TOOL_PATHS['star'], '--genomeDir', str(index_dir),
            '--runThreadN', str(num_threads),
            '--readFilesIn', str(processed_reads),
            '--outFileNamePrefix', str(output_dir / "EndToEnd_to_dc_"),
            '--outReadsUnmapped', 'Fastx',
            '--outSAMprimaryFlag', 'AllBestScore',
            '--alignEndsType', 'EndToEnd',
            '--scoreDelOpen', '0', '--scoreDelBase', '0',
            '--clip5pNbases', '10', '--alignIntronMin', '100'
        ]
        subprocess.run(cmd, check=True)
        return output_dir / "EndToEnd_to_dc_Aligned.out.sorted.bam"
        
    except subprocess.CalledProcessError as e:
        raise PipelineError(f"Alignment failed: {e.cmd}")

def find_gRNA_position(ref_sequence: str, gRNA_sequence: str, ref_coord: int) -> int:
    """Find gRNA position in reference sequence"""
    gRNA_rev_comp = str(Seq(gRNA_sequence).reverse_complement())
    pos_in_ref = ref_sequence.find(gRNA_rev_comp)
    
    if pos_in_ref == -1:
        raise ValueError("gRNA reverse complement not found in reference sequence")

    if pos_in_ref <= ref_coord < pos_in_ref + len(gRNA_rev_comp):
        pos_in_rev_gRNA = ref_coord - pos_in_ref
        return len(gRNA_rev_comp) - pos_in_rev_gRNA + 1
    return 'Unedited'

def get_switching_sites(cigar: str, refstart: int, edit_region_start: int, edit_region_end: int) -> List[Tuple]:
    """Extract deletions from CIGAR string"""
    cigar_pat = re.compile(r'\d+[A-Z]')
    tokens = cigar_pat.findall(cigar)
    refpos = refstart
    deletions = []
    
    for tok in tokens:
        op = tok[-1]
        length = int(tok[:-1])
        if op == 'D' and (edit_region_start <= refpos <= edit_region_end):
            deletions.append((refpos, refpos + length - 1, length))
        if op in ('M', 'D', 'N'):
            refpos += length
    return deletions

def analyze_deletions(
    bam_file: Path,
    reference_seq: str,
    guide_rna: str,
    sample_name: str,
    output_dir: Path
) -> None:
    """Integrated deletion analysis from original script2"""
    try:
        output_csv = output_dir / "df_level_version4.2.csv"
        gRNA_rc = str(Seq(guide_rna).reverse_complement())
        gRNA_start = reference_seq.find(gRNA_rc)
        gRNA_end = gRNA_start + len(gRNA_rc) - 1

        # Process BAM file
        max_value_list = []
        deletion_sites = defaultdict(int)
        no_deletion_count = 0
        read_count = 0

        with pysam.AlignmentFile(bam_file) as bam:
            for aln in tqdm(bam):
                if not aln.cigarstring:
                    continue
                
                read_count += 1
                deletions = get_switching_sites(
                    aln.cigarstring, aln.reference_start, 
                    gRNA_start, gRNA_end
                )
                
                if not deletions:
                    no_deletion_count += 1
                
                max_value_list.extend([length for _, _, length in deletions])
                for site in deletions:
                    deletion_sites[site] += 1

        # Write results
        with open(output_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['sample', 'start', 'end', 'length', 
                           'read_count', 'Total_read', 'del_freq','gRNA_start'])
            
            for (start, end, length), count in deletion_sites.items():
                del_freq = count / read_count
                grna_pos = find_gRNA_position(reference_seq, guide_rna, end+1)
                writer.writerow([sample_name, start, end, length, 
                               count, read_count, del_freq, grna_pos])
            
            # Write no-deletion entry
            del_freq = no_deletion_count / read_count
            writer.writerow([sample_name, 'no_deletion', 'no_deletion', 
                           'no_deletion', no_deletion_count, read_count, 
                           del_freq, 'Unedited'])

    except Exception as e:
        raise PipelineError(f"Deletion analysis failed: {str(e)}")




def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')



def main():
    """Main pipeline execution flow"""
    parser = argparse.ArgumentParser(
        description='CRISPR Sequence Analysis Pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--sample_name', required=True, help='Unique sample identifier')
    parser.add_argument('--read1', type=Path, required=True, help='Path to read 1 FASTQ')
    parser.add_argument('--read2', type=Path, required=True, help='Path to read 2 FASTQ')
    parser.add_argument('--raw_seq', required=True, help='Reference DNA sequence')
    parser.add_argument('--gRNA', required=True, help='Guide RNA sequence')
    parser.add_argument('--target_in_read1', type=str2bool, help='')
    
    # Optional arguments
    parser.add_argument('--work_dir', type=Path, default=Path.cwd(), help='Working directory')
    parser.add_argument('--threads', type=int, default=5, help='Number of processing threads')
    parser.add_argument('--steps', type=int, choices=[1,12,123,3,4], required= False,
                      help='Processing steps to execute (1=index, 2=process, 3=align, 4=analyze)')
    
    args = parser.parse_args()
    
    try:
        # Validate environment
        validate_tools(TOOL_PATHS)
        
        # Setup directory structure
        sample_dir = create_sample_directories(args.work_dir, args.sample_name)
        
        if args.target_in_read1:
            read_input = args.read1
        else:
            read_input = args.read2
        
        # Execute processing steps
        if args.steps is not None :
            print('Running')
            if args.steps in [12, 123]:
                processed_reads = process_reads(
                    read_input, sample_dir,
                )
            
            if args.steps in [1, 123]:
                build_star_index(args.raw_seq.upper(), sample_dir)
            
            if args.steps == 123:
                bam_file = run_star_alignment(processed_reads, sample_dir, sample_dir, args.threads)
                analyze_deletions(bam_file, args.raw_seq.upper(), 
                                args.gRNA.upper(), args.sample_name, sample_dir)
        else:
            print('Running all steps')
            processed_reads = process_reads(read_input, sample_dir,)
            build_star_index(args.raw_seq.upper(), sample_dir)
            bam_file = run_star_alignment(processed_reads, sample_dir, sample_dir, args.threads)
            analyze_deletions(bam_file, args.raw_seq.upper(), args.gRNA.upper(), args.sample_name, sample_dir)

    except PipelineError as e:
        print(f"Pipeline failed: {str(e)}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {str(e)}", file=sys.stderr)
        sys.exit(2)

if __name__ == "__main__":
    main()
