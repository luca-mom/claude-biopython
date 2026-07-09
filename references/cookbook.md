# Biopython Cookbook

Copy-paste-ready snippets, organized by task. Adapt variable names/paths to the user's actual files.

## Table of contents
1. Sequence I/O (FASTA/FASTQ/GenBank)
2. Sequence manipulation
3. Pairwise & multiple sequence alignment
4. BLAST (local and NCBI web)
5. Entrez (NCBI database queries)
6. Restriction enzyme analysis
7. Protein structures (PDB/mmCIF)
8. Phylogenetics

---

## 1. Sequence I/O

### Read FASTA
```python
from Bio import SeqIO

for record in SeqIO.parse("input.fasta", "fasta"):
    print(record.id, len(record.seq))
```

### Read FASTQ (with quality scores)
```python
from Bio import SeqIO

for record in SeqIO.parse("input.fastq", "fastq"):
    quals = record.letter_annotations["phred_quality"]
    print(record.id, len(record.seq), min(quals), max(quals))
```

### Load all records into a dict (small files only)
```python
records = SeqIO.to_dict(SeqIO.parse("input.fasta", "fasta"))
seq = records["some_id"].seq
```

### Write / convert between formats
```python
from Bio import SeqIO
SeqIO.convert("input.gb", "genbank", "output.fasta", "fasta")

# Or manually filter then write
records = (r for r in SeqIO.parse("input.fasta", "fasta") if len(r.seq) > 100)
SeqIO.write(records, "filtered.fasta", "fasta")
```

### Parse GenBank (includes annotations/features)
```python
for record in SeqIO.parse("input.gbk", "genbank"):
    print(record.id, record.description)
    for feature in record.features:
        if feature.type == "CDS":
            print(feature.location, feature.qualifiers.get("gene"))
```

---

## 2. Sequence manipulation

```python
from Bio.Seq import Seq

seq = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

seq.reverse_complement()          # DNA reverse complement
seq.transcribe()                  # DNA -> RNA (T -> U)
seq.translate()                   # -> protein (stops at first stop codon by default: table=1)
seq.translate(table=11, to_stop=True)   # bacterial code, stop at first stop, no '*' appended

from Bio.SeqUtils import gc_fraction
gc_fraction(seq) * 100            # %GC

from Bio.SeqUtils import molecular_weight
molecular_weight(seq, seq_type="DNA")
```

### Mutable sequences
```python
from Bio.Seq import MutableSeq
mseq = MutableSeq(str(seq))
mseq[0] = "T"
```

### Find ORFs / codons manually
```python
codon_table = {
    "start_codons": ["ATG"],
    "stop_codons": ["TAA", "TAG", "TGA"],
}
# For real ORF finding, prefer scanning all 6 frames (3 forward + 3 reverse-complement)
# and translating with to_stop=False to locate '*' positions, or use a dedicated tool
# like Prodigal/ORFfinder for genome-scale work — Biopython alone is not an ORF caller.
```

---

## 3. Pairwise & multiple sequence alignment

### Pairwise alignment (modern API — `Bio.Align.PairwiseAligner`)
```python
from Bio.Align import PairwiseAligner

aligner = PairwiseAligner()
aligner.mode = "global"          # or "local" for Smith-Waterman-style
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -2
aligner.extend_gap_score = -0.5

alignments = aligner.align("ACCGGT", "ACG")
best = alignments[0]
print(best)
print(best.score)
```

Note: `Bio.pairwise2` is deprecated in recent Biopython versions — use `Bio.Align.PairwiseAligner` instead.

### Reading a multiple sequence alignment file
```python
from Bio import AlignIO

alignment = AlignIO.read("input.aln", "clustal")
for record in alignment:
    print(record.id, record.seq)

print(alignment.get_alignment_length())
```

### Running an external MSA tool (Clustal Omega / MUSCLE)
Biopython does not compute MSAs itself for more than pairwise sequences — it wraps external tools:
```python
from Bio.Align.Applications import ClustalOmegaCommandline
# Requires clustalo installed on the system (not a pip package)
cmdline = ClustalOmegaCommandline(infile="input.fasta", outfile="aligned.fasta", verbose=True, auto=True)
stdout, stderr = cmdline()
```
Check the external binary is actually installed (`which clustalo`) before relying on this — it's a system dependency, not something `pip install biopython` provides.

---

## 4. BLAST

### NCBI web BLAST (no local database needed, but slow — network round trip)
```python
from Bio.Blast import NCBIWWW, NCBIXML

result_handle = NCBIWWW.qblast("blastn", "nt", open("query.fasta").read())

with open("blast_results.xml", "w") as out:
    out.write(result_handle.read())

# Parse results
with open("blast_results.xml") as f:
    blast_record = NCBIXML.read(f)

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < 0.01:
            print(alignment.title[:80])
            print("e-value:", hsp.expect, "identity:", hsp.identities, "/", hsp.align_length)
```
Warn the user: NCBI web BLAST can take minutes and NCBI rate-limits/blocks abusive automated use — this isn't suitable for looping over many queries. For bulk work, recommend local BLAST+.

### Local BLAST+ (requires `blastn`/`blastp` binaries installed separately, not via pip)
```python
from Bio.Blast.Applications import NcbiblastnCommandline

blastn_cline = NcbiblastnCommandline(
    query="query.fasta", db="my_local_db", evalue=0.001,
    outfmt=5, out="results.xml"
)
stdout, stderr = blastn_cline()

from Bio.Blast import NCBIXML
with open("results.xml") as f:
    record = NCBIXML.read(f)
```

---

## 5. Entrez (NCBI database queries)

**Always set `Entrez.email` first — required by NCBI.**

```python
from Bio import Entrez, SeqIO

Entrez.email = "your_email@example.com"  # tell user to replace with a real address

# Fetch a nucleotide sequence by accession
handle = Entrez.efetch(db="nucleotide", id="NM_001301717", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()
print(record.description, len(record.seq))

# Search for records matching a term
handle = Entrez.esearch(db="pubmed", term="CRISPR AND 2024[pdat]", retmax=20)
search_results = Entrez.read(handle)
handle.close()
id_list = search_results["IdList"]

# Fetch PubMed abstracts
handle = Entrez.efetch(db="pubmed", id=id_list, rettype="abstract", retmode="text")
print(handle.read())
handle.close()
```

If the user has an NCBI API key, set `Entrez.api_key = "..."` to raise the rate limit from 3 to 10 requests/second.

---

## 6. Restriction enzyme analysis

```python
from Bio.Seq import Seq
from Bio.Restriction import EcoRI, BamHI, AllEnzymes, Analysis

seq = Seq("GAATTCGGATCCAAGCTT")

EcoRI.search(seq)          # cut positions for a single enzyme

# Analyze with many enzymes at once
analysis = Analysis(AllEnzymes, seq)
result = analysis.with_sites()   # dict of enzyme -> cut positions, only enzymes that cut
```

---

## 7. Protein structures (PDB/mmCIF)

### Parse and inspect a structure
```python
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)   # QUIET suppresses non-fatal format warnings
structure = parser.get_structure("my_protein", "structure.pdb")

for model in structure:
    for chain in model:
        print("Chain:", chain.id, "residues:", len(list(chain)))
        for residue in chain:
            if residue.id[0] == " ":  # skip heteroatoms/waters
                pass
```

### mmCIF (preferred for modern PDB entries)
```python
from Bio.PDB import MMCIFParser
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("my_protein", "structure.cif")
```

### Download a structure directly from the PDB
```python
from Bio.PDB import PDBList
pdbl = PDBList()
pdbl.retrieve_pdb_file("1CRN", pdir=".", file_format="pdb")
```

### Measure distance between atoms / basic geometry
```python
atom1 = structure[0]["A"][10]["CA"]
atom2 = structure[0]["A"][50]["CA"]
distance = atom1 - atom2   # Biopython overloads subtraction for Euclidean distance
```

### Compute secondary structure / solvent accessibility
Requires the external `DSSP` binary installed separately:
```python
from Bio.PDB import DSSP
model = structure[0]
dssp = DSSP(model, "structure.pdb")  # needs `mkdssp` on PATH
```

---

## 8. Phylogenetics

### Read and draw a tree
```python
from Bio import Phylo

tree = Phylo.read("tree.nwk", "newick")
Phylo.draw_ascii(tree)              # quick terminal view

# Matplotlib rendering (for a saved image)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 10))
Phylo.draw(tree, axes=ax, do_show=False)
fig.savefig("tree.png")
```

### Build a tree from a distance matrix (simple neighbor-joining / UPGMA)
```python
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

alignment = AlignIO.read("aligned.fasta", "fasta")
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)

constructor = DistanceTreeConstructor()
nj_tree = constructor.nj(dm)     # neighbor-joining
upgma_tree = constructor.upgma(dm)

Phylo.write(nj_tree, "nj_tree.nwk", "newick")
```
