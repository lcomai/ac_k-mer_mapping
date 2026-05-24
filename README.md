# Aho-Corasick k-mer genome mapper

Map exact instances of a list of k-mers onto a reference genome FASTA using the
Aho-Corasick algorithm. The mapper reports every matching position on the
forward strand and, by default, the reverse-complement strand.

This is a simple utility developed for mapping exact instances of k-mers on
small to medium plant genomes. It is intended for straightforward exploratory
and reproducible analyses, not as a full-featured read aligner or general genome
indexing system.

The program does not require chromosome-style names. FASTA records can be named
`Chr1`, `1`, `NC_003070.9`, `contig_42`, or anything else. The record name is
written to the output as `seq_id`; it is not interpreted during mapping.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

## Quick start

```bash
python3 map_kmers_ac.py \
  /Users/lucacomai/Documents/Data_analysis_23/kmers/busco_kmers.txt \
  /Users/lucacomai/Documents/Data_analysis_23/Genomes/Col-PEK1.5/file1.Col-PEK1.5_Chr1-5_20220523.fasta \
  outputs/busco_kmers_col_pek15.tsv
```

The genome argument can be either a single FASTA file or a directory containing
FASTA files with suffixes `.fa`, `.fasta`, `.fna`, `.fas`, or `.fsa`.

If a directory contains more than one FASTA, the mapper stops and asks you to be
explicit. This avoids accidentally mapping against multiple assemblies, such as
an old and current genome stored in the same directory. Use an exact FASTA path,
`--fasta-glob`, or `--all-fastas`.

Outputs:

- `outputs/busco_kmers_col_pek15.tsv`: all mapped k-mer instances
- `outputs/busco_kmers_col_pek15.counts.tsv`: optional counts per k-mer
- `outputs/busco_kmers_col_pek15.tsv.no_match`: valid k-mers with no hit
- `outputs/busco_kmers_col_pek15.tsv.invalid.tsv`: invalid input k-mer lines
- `outputs/busco_kmers_col_pek15.tsv.summary.tsv`: run statistics

Output columns:

| column | meaning |
| --- | --- |
| `kmer` | original input k-mer |
| `seq_id` | FASTA record name |
| `start` | 1-based inclusive start coordinate |
| `end` | 1-based inclusive end coordinate |
| `strand` | `+` for input orientation, `-` for reverse-complement hit |
| `source_fasta` | FASTA file that contained the hit |

## Optional step: count mapped positions per k-mer

After mapping, tally the number of genome positions found for each k-mer:

```bash
python3 tally_kmer_hits.py \
  outputs/busco_kmers_col_pek15.tsv \
  outputs/busco_kmers_col_pek15.counts.tsv
```

The count output has one row per mapped k-mer:

| column | meaning |
| --- | --- |
| `kmer` | original input k-mer |
| `mapped_positions` | number of mapped instances in the genome |

To keep the original input order and include k-mers with no mapped positions,
provide the original k-mer file:

```bash
python3 tally_kmer_hits.py \
  outputs/busco_kmers_col_pek15.tsv \
  outputs/busco_kmers_col_pek15.counts.tsv \
  --kmers /Users/lucacomai/Documents/Data_analysis_23/kmers/busco_kmers.txt
```

## FASTA record names

By default, the mapper uses the first token after `>` as `seq_id`, matching
common FASTA behavior:

```text
>Chr1 Col-PEK1.5 chromosome 1
```

is reported as `Chr1`.

If you need the full header in the output, use:

```bash
python3 map_kmers_ac.py kmers.txt reference.fa out.tsv --name-mode full
```

This makes the script robust to different genome naming styles because mapping
is performed only on sequence content. Name differences only affect labels in
the output and any downstream program that expects a specific name such as
`Chr1` instead of `1`.

## Input k-mers

The input file is plain text. The first whitespace-delimited field on each
non-empty, non-comment line is read as the k-mer.

```text
ACGTTGCAACGTTGCAACGTTGC
TTGCAACGTTGCAACGTTGCAAC
```

K-mers can have mixed lengths. By default, only `A`, `C`, `G`, and `T` are
accepted. Use `--allow-n` to allow `N` as a literal character. `N` is not treated
as a wildcard.

Very short k-mers can occur millions of times in a genome. For production runs,
use biologically meaningful k-mer lengths, such as the 23-mers in the BUSCO
k-mer list.

## Useful options

```bash
# Scan FASTA files inside nested genome directories
python3 map_kmers_ac.py kmers.txt genome_dir out.tsv --recursive --all-fastas

# Select FASTAs from a mixed directory
python3 map_kmers_ac.py kmers.txt genome_dir out.tsv --fasta-glob '*Col-PEK1.5*.fasta'

# Search only input orientation, not reverse complements
python3 map_kmers_ac.py kmers.txt reference.fa out.tsv --forward-only

# Preserve full FASTA headers
python3 map_kmers_ac.py kmers.txt reference.fa out.tsv --name-mode full
```

## Development check

```bash
python3 -m unittest discover tests
```

## Notes for GitHub release

Choose and add a license before publishing if other people should be able to
reuse or modify the code under explicit terms.
