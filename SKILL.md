---
name: biopython
description: Use this skill for any bioinformatics task in Python involving DNA/RNA/protein sequences, FASTA/FASTQ/GenBank/PDB files, sequence alignment, BLAST searches, motif finding, phylogenetics, transcription/translation, restriction analysis, or protein structure parsing. Trigger whenever the user mentions Biopython, .fasta/.fa/.fastq/.gb/.gbk/.pdb/.cif files, sequence I/O, BLAST, multiple sequence alignment, Entrez/NCBI lookups, codon translation, GC content, or working with biological sequence/structure data in Python — even if they just say "parse this sequence file" or "analyze this protein structure" without naming the library.
---

# Biopython

Biopython is the standard Python toolkit for computational biology: parsing sequence/structure file formats, manipulating DNA/RNA/protein sequences, running and parsing BLAST, working with alignments and phylogenetic trees, and querying NCBI databases.

## Setup

Check if installed before writing code:
```bash
python3 -c "import Bio; print(Bio.__version__)" || pip install biopython --break-system-packages
```

## Core workflow

1. **Identify the file format / task type** from the user's request or the file extension (see table below).
2. **Read `references/cookbook.md`** for the relevant section — it has copy-paste-ready code for each common task (I/O, sequence ops, alignment, BLAST, Entrez, structures, phylogenetics).
3. **Write a script**, not one-off interpreter commands, for anything beyond a single trivial lookup — bioinformatics scripts are reused and inspected.
4. **Never guess at biological results.** If asked for a fact like a gene's function, organism source, or database annotation, fetch it via `Bio.Entrez` or say you don't have it — don't fabricate accession numbers, sequences, or annotations.
5. **Large files** (>~50MB FASTQ/BAM-scale data): use `SeqIO.parse` (a generator/iterator) rather than `SeqIO.to_dict`, to avoid loading everything into memory. Never call `list()` on a parse iterator for large files.

## Format quick-reference

| Extension | Format | Bio module |
|---|---|---|
| .fasta, .fa, .fna, .faa | FASTA | `Bio.SeqIO` |
| .fastq, .fq | FASTQ (with quality scores) | `Bio.SeqIO` |
| .gb, .gbk, .genbank | GenBank | `Bio.SeqIO` |
| .embl | EMBL | `Bio.SeqIO` |
| .aln, .clustal | Clustal alignment | `Bio.AlignIO` |
| .sto, .stockholm | Stockholm alignment | `Bio.AlignIO` |
| .phy, .phylip | Phylip alignment | `Bio.AlignIO` |
| .nwk, .nexus, .nex | Phylogenetic tree | `Bio.Phylo` |
| .pdb | Protein structure (legacy) | `Bio.PDB.PDBParser` |
| .cif, .mmcif | Protein structure (modern) | `Bio.PDB.MMCIFParser` |
| BLAST XML output | BLAST results | `Bio.Blast.NCBIXML` |

## Key gotchas (read before writing code)

- **`Seq` objects are mostly immutable** in modern Biopython. To edit in place, convert with `MutableSeq` or use string operations and re-wrap in `Seq(...)`.
- **Always set `Entrez.email`** before any `Bio.Entrez` call — NCBI requires it and may block requests without it. Use a placeholder like `"your_email@example.com"` and tell the user to replace it with their real address; NCBI's terms require a genuine contact email for API use.
- **`SeqIO.parse` returns an iterator**, consumed once. If you need to loop over records twice, either re-parse the file or convert to a list first (only for small files).
- **Translation requires the correct reading frame and table.** Default `.translate()` uses the standard genetic code (table 1); mitochondrial/bacterial sequences need `table=2`, `table=11`, etc. Ask the user if the organism/context is ambiguous.
- **PDB parsing is noisy by default** — wrap `PDBParser(QUIET=True)` unless the user wants to see structure warnings.
- **Coordinates in Biopython are 0-indexed, half-open**, matching Python slicing — this differs from 1-indexed conventions used in some bio file formats (e.g., GFF, VCF are 1-indexed). Be explicit about which convention is in play when converting between formats.

## When to consult the cookbook

Open `references/cookbook.md` for ready-to-adapt code whenever the task involves:
- Reading/writing/converting sequence files (FASTA, FASTQ, GenBank)
- Sequence manipulation (reverse complement, transcription, translation, GC content)
- Running BLAST (local or NCBI web) and parsing results
- Pairwise or multiple sequence alignment
- Fetching data from NCBI (Entrez: nucleotide, protein, PubMed)
- Restriction enzyme analysis
- Parsing/analyzing protein structures (PDB/mmCIF)
- Building or rendering phylogenetic trees

Don't reproduce the whole cookbook file to the user — read what's relevant, adapt it to their actual data, and present a working script or direct answer.
