# 🧬 bioseqkit
*A lightweight Python toolkit for DNA/RNA operations and FASTQ read filtering,
completed by instruments for bioinformatics files processing*

---

## 📘 Overview
**bioseqkit** is a simple educational package that implements essential bioinformatics utilities:
- DNA/RNA sequence manipulations: transcription, reverse, complement, reverse-complement.  
- FASTQ-like read filtering based on GC content, sequence length, and Phred + 33 quality.
- FASTQ, FASTA and GBK files manipulation
  

All functionality is implemented **without external dependencies** — only the Python 3.10+ standard library.

---

## 📁 Repository structure
```
bioseqkit/
├─ README.md
├─ bioseqkit.py                # main entry point (imports + 2 functions)
├─ bio_files_processor.py      # additional file with instruments for files processing
└─ modules/                    # helper modules
   ├─ dna_rna.py               # DNA/RNA sequence utilities
   ├─ fastq.py                 # FASTQ filtering helpers
   ├─ bio_files.py             # Files processing helpers
   └─ types_alliases.py        # shared type definitions
```

---

## ⚙️ Installation 
Clone the repository and navigate into it:
```bash
git clone https://github.com/arnautoleg/bioseqkit.git
cd bioseqkit
git switch hw5-Files # before the pull request
```

No installation of extra packages is required.

---

## 🚀 Quick start

### 1️⃣ DNA/RNA utilities
Use the function `run_dna_rna_tools(*args)` to apply one of the supported procedures:

```python
from bioseqkit import run_dna_rna_tools

run_dna_rna_tools('TTUU', 'is_nucleic_acid') # False !!
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATc', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']

```

**Supported operations:**
- `"is_nucleic_acid"`
- `"transcribe"`
- `"reverse"`
- `"complement"`
- `"reverse_complement"`

---

### 2️⃣ FASTQ filtering
Use `filter_fastq(seqs, gc_bounds, length_bounds, quality_threshold)`  
to keep only reads that satisfy all filters.

```python
from bioseqkit import filter_fastq

out = filter_fastq(
    input_fastq="example_data/example_fastq.fastq",
    output_fastq="example_data/example_filtered.fastq",
    gc_bounds=(40, 60),
    length_bounds=(50, 300),
    quality_threshold=20
)

```
All filtrated reads are saved in "example_data/example_filtered.fastq" directory

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

## 🧩 Module summary

| Module | Description |
|---------|--------------|
| `dna_rna.py` | Sequence validation, transcription, reverse/complement operations |
| `fastq.py` | Filtering helpers: GC%, length, quality, read_fastq_to_dict, write_fastq_from_dict, filter_fastq_dict |
| `bio_files.py` | Files processing helpers: _matches_anchor, _feature_header, _wrap60, _flush_current_feature, _finish_record |
| `types_alliases.py` | Common type aliases (`seq_tuple`, `seq_dict`, `bounds`) |

---

## 🧪 Testing (in development)
You can add or run tests via [pytest](https://pytest.org):
```bash
pytest -q
```
Recommended minimal test set includes checks for:
- GC % and length inclusion.
- Upper-bound interpretation.
- Phred + 33 quality threshold behavior.
- Correct reverse/complement output.

---

## 📄 License
This project is distributed under the **MIT License**.

---

## ✉️ Author
Developed by **Oleg Arnaut**,  
Bioinformatic Institute
