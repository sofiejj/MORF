# MORF

**MORF: A Computational Framework for Synonymous Optimisation of Co-Encoding Regions in Protein-Coding Sequences**

MORF (**multiple overlapping reading frames**) is a computational framework for identifying and optimising candidate co-encoding regions (CERs) between protein-coding nucleotide sequences.

The framework uses synonymous codon substitutions to increase nucleotide-level agreement between overlapping coding sequences while preserving the amino-acid sequence encoded by each input sequence.

This repository contains the MATLAB implementation used for the proof-of-concept evaluation described in the accompanying preprint.

## Overview

MORF consists of three main stages:

1. **Candidate CER identification**  
   Sliding Hamming similarity (SHS) is evaluated across admissible relative nucleotide offsets between two protein-coding sequences.

2. **Synonymous codon optimisation**  
   The highest-ranking candidate CERs are independently optimised using a sequential two-pass procedure. In the first pass, sequence B is modified while sequence A is held fixed. The resulting sequence B is then held fixed while sequence A is modified in the second pass. Synonymous substitutions are retained only where they produce a strict increase in local nucleotide-level agreement.

3. **Final CER selection**  
   Candidate regions are compared according to their post-optimisation nucleotide agreement, and the highest-scoring configuration is selected.

Only synonymous codon substitutions are permitted in the current implementation. The amino-acid sequences encoded by both input sequences are therefore preserved.

## Proof-of-Concept Sequences

The evaluation presented in the accompanying preprint uses Green Fluorescent Protein (GFP) and DsRed.

| Protein | GenBank accession | Coding-sequence length |
|---------|--------------------|------------------------|
| GFP     | M62653.1           | 717 nt                 |
| DsRed   | AF168419.2         | 678 nt                 |

The exact nucleotide sequences used for the evaluation are contained in `seq_chooser.m`.

## Implementation

The current MORF prototype is implemented in MATLAB.

For the GFP–DsRed proof-of-concept evaluation, the following settings are specified in `mainB.m`:

- Sequence A: GFP
- Sequence B: DsRed
- Minimum CER length: 339 nt
- Number of candidate CERs retained: 10
- Reading-frame configurations considered: forward RF 0, RF 1 and RF 2
- Permitted substitutions: synonymous codons only
- Optimisation order: sequence B followed by sequence A
- Acceptance criterion: strict increase in nucleotide agreement

The optimisation is greedy and sequential. It does not perform exhaustive enumeration of all possible synonymous sequence combinations and therefore does not guarantee identification of the globally optimal synonymous representation.

## Running MORF

Clone or download this repository and set `MORF_v1` as the current MATLAB working directory.

Run:

```matlab
mainB
```

`mainB.m` performs the complete GFP–DsRed proof-of-concept analysis using the sequences and parameters described above.

The implementation:

- evaluates sliding Hamming similarity across admissible relative offsets;
- retains the ten highest-scoring candidate CERs;
- performs sequential two-pass synonymous codon optimisation on each candidate;
- calculates nucleotide agreement before and after optimisation;
- reinserts each optimised CER into the corresponding full-length coding sequences;
- translates the resulting sequences to verify amino-acid sequence preservation; and
- generates diagnostic and summary visualisations.

## Reproducing the Preprint Results

Running `mainB.m` with the supplied GFP and DsRed sequences and default settings reproduces the proof-of-concept analysis reported in the accompanying preprint.

Across the ten retained candidate CERs, mean nucleotide agreement increases from **32.06%** before optimisation to **53.43%** following optimisation.

The highest-scoring post-optimisation candidate is **CER 3**, for which nucleotide agreement increases from **32.17% to 60.05%**.

**CER 7** exhibits the largest absolute improvement, increasing from **31.95% to 60.00%**, corresponding to an increase of **28.05 percentage points**.

Following optimisation, the translated GFP and DsRed amino-acid sequences are identical to their corresponding original sequences for all ten candidate CERs.

## Repository Structure

```text
MORF_v1/
├── mainB.m
├── seq_chooser.m
├── slidingHamming.m
├── build_numeric_codon_table.m
├── find_best_codon_pair.m
├── map_base2num.m
├── map_num2base.m
└── translate_nt.m
```

### File descriptions

- **`mainB.m`** — main MORF workflow and entry point for reproducing the proof-of-concept analysis.
- **`seq_chooser.m`** — contains the GFP and DsRed nucleotide sequences used in the evaluation.
- **`slidingHamming.m`** — evaluates nucleotide-level sliding Hamming similarity and identifies the leading candidate CERs.
- **`build_numeric_codon_table.m`** — constructs the standard genetic-code lookup used to identify synonymous codon alternatives.
- **`find_best_codon_pair.m`** — evaluates candidate codon alternatives according to local nucleotide agreement.
- **`map_base2num.m`** — converts nucleotide characters to the numeric representation used internally by MORF.
- **`map_num2base.m`** — converts the internal numeric representation back to nucleotide characters.
- **`translate_nt.m`** — translates nucleotide sequences using the standard genetic code for verification of amino-acid sequence preservation.

## Requirements

- MATLAB

The current implementation uses standard MATLAB functionality and does not require input data files: the GFP and DsRed sequences used for the proof-of-concept analysis are included directly in `seq_chooser.m`.

## Current Scope

MORF is currently a proof-of-concept computational implementation.

The present version:

- considers pairwise co-encoding of two protein-coding sequences;
- considers forward-orientation reading-frame configurations;
- restricts sequence modification to synonymous codon substitutions; and
- evaluates candidate CERs primarily according to nucleotide-level agreement.

The current results demonstrate computational sequence compatibility and amino-acid sequence preservation. They do not by themselves establish successful biological expression or function of an optimised co-encoding construct.

Future development will include evaluation across additional protein pairs, reverse-complement and antisense configurations, additional biologically relevant optimisation constraints, and experimental validation.

## Citation

If you use MORF in your work, please cite:

> van de Bijlmer, S. (2026). *MORF: A Computational Framework for Synonymous Optimisation of Co-Encoding Regions in Protein-Coding Sequences*. Preprint.

A DOI and preprint link will be added here following publication of the preprint.

## Author

**Sofie van de Bijlmer**

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.
