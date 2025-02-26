# SCISSOR: Selective Cleavages and Intramolecular Stitches of RNA

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)

This repository contains the analysis pipeline for the paper:  
**"Type III CRISPR-mediated flexible RNA excision with engineered guide RNAs"**  
*Molecular Cell (2025)*

## 📖 Overview
SCISSOR (Selective Cleavages and Intramolecular Stitches of RNA) is a novel RNA editing tool that enables precise RNA fragment excision using engineered guide RNAs with bulge loops. This repository contains the bioinformatics pipeline for analyzing CRISPR sequencing data and characterizing RNA deletion patterns.

Key features:
- Adapter trimming and sequence filtering
- Custom STAR index construction
- End-to-end sequence alignment
- Deletion pattern analysis
- Frameshift mutation characterization

## 📦 Dependencies
### Core Requirements
- Python 3.8+
- Cutadapt (v3.7)
- STAR (v2.7.10a)

## 🚀 Quick Start

pleace check the configuration for external tool paths first
1. Paired-end FASTQ files
2. Reference sequence (DNA)
3. Guide RNA sequence without bulge loop

### Basic Usage
```bash
# After confirming the paths of each software in the python file
python scissor_call_deletion.py \
      --sample {sample} \
      --read1 {Path}/{sample}/{sample}_raw_1.fq.gz \
      --read2 {Path}/{sample}/{sample}_raw_2.fq.gz \
      --raw_seq {target_seq}\
      --gRNA {gRNA_seq} \
      --threads 10 \
      --target_in_read1 {if_target_in_read1} "
```
### Command Line Options
| Parameter          | Required | Description                          |
|---------------------|----------|--------------------------------------|
| `--sample_name`     | Yes      | Unique experiment identifier        |
| `--read1`           | Yes      | Path to Read 1 FASTQ file           |
| `--read2`           | Yes      | Path to Read 2 FASTQ file           |
| `--raw_seq`         | Yes      | Reference DNA sequence    |
| `--gRNA`            | Yes      | Guide RNA sequence    |
| `--target_in_read1` | Yes      | Target strand in Read1 (True/False) |
| `--work_dir`        | No       | Output directory (default: current) |
| `--threads`         | No       | CPU threads (default: 5)           |
| `--steps`           | No       | Pipeline steps (1=index, 2=process, 3=align, 4=analyze) |




### Pipeline Steps
1. **Adapter Processing**  
   - Remove specified sequences (`--sequence_remove`)
   - Retain target sequences (`--sequence_retain`)
   - Quality trimming and adapter removal

2. **Index Construction**  
   Build custom STAR index from reference sequence

3. **Sequence Alignment**  
   End-to-end mapping with optimized parameters

4. **Deletion Analysis**  
   Characterize deletion patterns and frameshift mutations



   

## 📄 Output Files
| File | Description |
|------|-------------|
| `reads.cutadapt.fastq` | Processed sequencing reads |
| `EndToEnd_to_dc_Aligned.out.bam` | Alignment results |
| `df_level_version4.2.csv` | Deletion analysis report |
| `std_out_*` | Process logs |

## 📚 Citation
```bibtex
Sun, Y. et al. Type III CRISPR-mediated flexible RNA excision with engineered guide RNAs. Mol. Cell (2025) doi:10.1016/j.molcel.2025.01.021.
```


## 📧 Contact
For technical support or scientific collaboration:  
Rui Zhang - zhangrui3@mail.sysu.edu.cn  
Yuanfan Sun - sunyf36@mail2.sysu.edu.cn

---

*Developed by Zhang Lab @ Sun Yat-Sen University*  
*Guangzhou, China | 2025*












