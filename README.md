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
- seqkit (v2.3.0)
- STAR (v2.7.10a)

## 🚀 Quick Start
### Input Files
1. Paired-end FASTQ files
2. Reference sequence (DNA)
3. Guide RNA sequence

### Basic Usage
```bash
python scissor_pipeline.py \
  --sample_name TEST01 \
  --read1 sample_R1.fastq \
  --read2 sample_R2.fastq \
  --raw_seq "ATGCGTA...TACGCTA" \
  --gRNA "GAGCTCG...AAACGGT" \
  --work_dir ./analysis \
  --threads 8
```

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
