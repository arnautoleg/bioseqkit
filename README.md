# 🧬 bioseqkit
*A lightweight Python toolkit for DNA/RNA operations and FASTQ read filtering.*

---

## 📘 Overview
**bioseqkit** is a simple educational package that implements essential bioinformatics utilities:
- DNA/RNA sequence manipulations: transcription, reverse, complement, reverse-complement.  
- Validation of nucleic acid sequences.  
- FASTQ-like read filtering based on GC content, sequence length, and Phred + 33 quality.  

All functionality is implemented **without external dependencies** — only the Python 3.10+ standard library.

---

## 📁 Repository structure
```
bioseqkit/
├─ README.md
├─ bioseqkit.py          # main entry point (imports + 2 functions)
└─ modules/              # helper modules
   ├─ dna_rna.py         # DNA/RNA sequence utilities
   ├─ fastq.py           # FASTQ filtering helpers
   └─ types_alliases.py  # shared type definitions
```

---

## ⚙️ Installation
Clone the repository and navigate into it:
```bash
git clone https://github.com/arnautoleg/bioseqkit.git
cd bioseqkit
```

Run directly with Python 3.10 or higher:
```bash
python bioseqkit.py
```

No installation of extra packages is required.

---

## 🚀 Quick start

### 1️⃣ DNA/RNA utilities
Use the function `run_dna_rna_tools(*args)` to apply one of the supported procedures:

```python
from bioseqkit import run_dna_rna_tools

run_dna_rna_tools("ATGC", "reverse")                 # -> 'CGTA'
run_dna_rna_tools("ATGC", "transcribe")              # -> 'AUGC'
run_dna_rna_tools("ATGC", "AUGC", "is_nucleic_acid") # -> [True, True]
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

seqs = {
    "read1": ("ATGCATGC", "IIIIIIII"),  # GC=50%, len=8, meanQ≈40
    "read2": ("ATATATAT", "IIIIIIII"),  # GC=0%,  len=8, meanQ≈40
    "read3": ("GCGCGC",   "!!!!!!"),    # GC=100%, len=6, meanQ≈0
}

filtered = filter_fastq(
    seqs,
    gc_bounds=(40, 100),     # GC % ∈ [40, 100]
    length_bounds=(6, 10),   # len ∈ [6, 10]
    quality_threshold=30,    # mean Q ≥ 30
)

print(filtered.keys())  # -> dict_keys(['read1'])
```

**Notes:**
- Passing a single number for `gc_bounds` or `length_bounds` treats it as the *upper bound* (0, x).  
- Bounds are *inclusive* on both sides.

---

## 🧩 Module summary

| Module | Description |
|---------|--------------|
| `dna_rna.py` | Sequence validation, transcription, reverse/complement operations |
| `fastq.py` | Filtering helpers: GC%, length, quality |
| `types_alliases.py` | Common type aliases (`SeqTuple`, `SeqDict`, `Bounds`) |

---

## 🧪 Testing
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
