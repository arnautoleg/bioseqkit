# 🧬 bioseqkit
*Educational project implementing biological sequence processing using
Object-Oriented Programming (OOP) and Biopython.*

---

## 📘 Overview
**bioseqkit** is a simple educational package that implements essential bioinformatics utilities:
- DNA/RNA sequence manipulations: transcription, reverse, complement, reverse-complement.  
- FASTQ-like read filtering based on GC content, sequence length, and Phred + 33 quality.
- FASTQ, FASTA and GBK files manipulation
  

External dependency: Biopython (see requirements.txt).

---

## 📁 Repository structure
```

bioseqkit/
├─ README.md
├─ bioseqkit.py                # OOP implementation of biological sequences
├─ bio_files_processor.py      # file processing utilities (FASTA, BLAST, GBK)
├─ cli.py                      # command-line interface (argparse)
├─ bioseqkit_test.py           # pytest tests
├─ requirements.txt            # project dependencies
├─ .gitignore                  # excludes venv, cache
└─ example_data/               # example input files

```

---

## ⚙️ Installation 
Clone the repository and navigate into it:
```bash
git clone https://github.com/arnautoleg/bioseqkit.git
cd bioseqkit
git switch HW21 # before the pull request
```


```bash
# Create and activate a virtual environment
python3 -m venv venv
source venv/bin/activate

# Install dependencies (IMPORTANT: use python -m pip)
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

# Check that Biopython is installed correctly
python -c "import Bio; print(Bio.__version__)"
```

---

## 🚀 Quick start

### 1️⃣ DNA/RNA utilities

All sequence operations are now implemented as class methods (DNASequence, RNASequence, AminoAcidSequence)

```python

from bioseqkit import DNASequence, is_nucleic_acid

# Transcription
print(DNASequence("ATG").transcribe().sequence)# AUG

# Reverse
print(DNASequence("ATG").reverse().sequence)# GTA

# Complement
print(DNASequence("AtG").complement().sequence)# TAC

# Reverse complement
print(DNASequence("ATc").reverse_complement().sequence)# GAT

# Validation
print(is_nucleic_acid("TTUU"))# False

```

**Supported sequence classes:**
- `"DNASequence"`
- `"RNASequence"`
- `"AminoAcidSequence"`

---

### 2️⃣ FASTQ filtering (Biopython)

Filtering is implemented with Biopython SeqIO / SeqRecord and supports 
read length, mean Phred quality, and GC%

```python
from bioseqkit import filter_fastq

out = filter_fastq(
    input_fastq="example_data/example_fastq.fastq",
    output_fastq="example_filtered.fastq",
    gc_bounds=(40, 60),        # GC% range
    length_bounds=(50, 300),   # length range
    quality_threshold=20,      # mean Phred (Phred+33)
    out_dir="filtered",        # output folder (created automatically)
)

```
Saved to filtered/example_filtered.fastq (by default).

###  3️⃣ Convert multi-line FASTA to one-line per record

```python

from bio_files_processor import convert_multiline_fasta_to_oneline

out = convert_multiline_fasta_to_oneline("example_data/example_multiline_fasta.fasta")
print("Written to:", out)  

```

###  4️⃣ Parse BLAST TXT and keep only best hits

```python

from bio_files_processor import parse_blast_best

hits = parse_blast_best(
    "example_data/example_blast_results.txt",
    "example_data/best_hits.txt",
)

```

###  5️⃣ Extract CDS neighbors from GenBank to FASTA

```python

from bio_files_processor import select_genes_from_gbk_to_fasta

select_genes_from_gbk_to_fasta(
    input_gbk="example_data/example_gbk.gbk",
    genes=["gltA", "sdhA", "sucA"],      
    output_fasta="example_data/flanks_example.faa",
    n_before=2,
    n_after=3,
)

```


---

## 🧪 Testing
You can add or run tests via [pytest](https://pytest.org):
```bash
pytest -q bioseqkit_test.py
```

The tests checks for:

1. DNA transcription
2. reverse and reverse-complement operations
3. amino acid residue fraction calculation
4. invalid alphabet handling (ValueError)
5. nucleic acid validation
6. bounds normalization
7. FASTA file conversion (input/output file test)
8. FASTQ filtering by quality
9. BLAST best-hit parsing
10. empty gene list validation

---

## 🖥️ Command Line Interface (CLI)

The tool can be used directly from the command line via `cli.py`.

```bash
python cli.py <command> [arguments]
```

- Use relative paths if you run commands from the project root
- Use absolute paths if files are located elsewhere

### Example: CLI help
You can inspect available arguments for each command using `--help`.

Example for FASTQ filtering:

```bash
python cli.py filter-fastq --help
```

 <img width="669" height="386" alt="image" src="https://github.com/user-attachments/assets/13de9090-5ac6-4208-b9a1-24ceb4d6b3cc" />



### Example convert multiline FASTA (absolute path in WSL):

```bash
python cli.py fasta-oneline \
  --input-fasta /mnt/c/Users/Admin/.../example_data/example_multiline_fasta.fasta \
  --output-fasta output.fasta
```
---

## 📝 Logging

The CLI logs execution details into a file: bioseqkit.log


<img width="1120" height="383" alt="image" src="https://github.com/user-attachments/assets/d5613140-4e53-447e-bfe5-d9ab9853ab9b" />



## 📄 License
This project is distributed under the **MIT License**.

---

## ✉️ Author
Developed by **Oleg Arnaut**,  
Bioinformatic Institute
