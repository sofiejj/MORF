# MORF

**MORF: A Computational Framework for Synonymous Optimisation of Co-Encoding Regions in Protein-Coding Sequences**

MORF (multiple overlapping reading frames) is a computational framework for
identifying and optimising candidate co-encoding regions (CERs) between
protein-coding nucleotide sequences.

The framework uses synonymous codon substitutions to increase nucleotide-level
agreement between overlapping coding sequences while preserving the amino-acid
sequence encoded by each input sequence.

This repository contains the MATLAB implementation used for the proof-of-concept
evaluation described in the accompanying preprint.

## Overview

MORF consists of three main stages:

1. **Candidate CER identification**  
   Sliding Hamming similarity (SHS) is evaluated across admissible relative
   nucleotide offsets between two protein-coding sequences.

2. **Synonymous codon optimisation**  
   The highest-ranking candidate CERs are independently optimised using a
   sequential two-pass procedure. Synonymous substitutions are accepted only
   where they produce a strict increase in local nucleotide-level agreement.

3. **Final CER selection**  
   Candidate regions are compared according to their post-optimisation
   nucleotide agreement, and the highest-scoring configuration is selected.

Only synonymous codon substitutions are permitted in the current implementation.
The native amino-acid sequences of both input proteins are therefore preserved.

## Proof-of-Concept Sequences

The evaluation presented in the accompanying preprint uses Green Fluorescent
Protein (GFP) and DsRed.

| Protein | GenBank accession | Coding-sequence length |
|---------|--------------------|------------------------|
| GFP     | M62653.1           | 717 nt                 |
| DsRed   | AF168419.2         | 678 nt                 |

The exact nucleotide sequences used for the evaluation are included in this
repository.

## Implementation

The current MORF prototype is implemented in MATLAB.

For the GFP-DsRed proof-of-concept evaluation, the following settings were used:

- Minimum CER length: 339 nt
- Number of candidate CERs retained: 10
- Reading-frame configurations considered: forward RF 0, RF 1 and RF 2
- Permitted substitutions: synonymous codons only
- Optimisation order: sequence B followed by sequence A
- Acceptance criterion: strict increase in nucleotide agreement

The optimisation is greedy and sequential. It does not perform exhaustive
enumeration of all possible synonymous sequence combinations and therefore does
not guarantee identification of the globally optimal synonymous representation.

## Running MORF

run('mainB.m')

The implementation identifies candidate CERs, performs synonymous codon
optimisation, evaluates nucleotide agreement before and after optimisation,
and generates diagnostic visualisations of the resulting candidate regions.

## Reproducing the Preprint Results

Using the GFP and DsRed sequences and settings described above, the implementation
retains ten candidate CERs for optimisation.

Across these candidates, mean nucleotide agreement increases from 32.06%
before optimisation to 53.43% following optimisation.

The highest-scoring post-optimisation candidate is CER 3, for which
nucleotide agreement increases from 32.17% to 60.05%.

The amino-acid sequences encoded by GFP and DsRed remain unchanged following
optimisation for all ten candidate CERs.

## Repository Structure

MORF_v1/
├── build_numeric_codon_table.m
├── find_best_codon_pair.m
├── mainB.m
├── map_base2num.m
├── map_num2base.m
├── seq_chooser.m
├── slidingHamming.m
└── translate_nt.m

## Requirements

MATLAB

## Citation

If you use MORF in your work, please cite:

van de Bijlmer, S. (2026). MORF: A Computational Framework for Synonymous
Optimisation of Co-Encoding Regions in Protein-Coding Sequences. Preprint.

A DOI and preprint link will be added here following publication of the preprint.

## Current Scope

MORF is currently a proof-of-concept computational implementation.

The present version:

- considers pairwise co-encoding of two protein-coding sequences;
- considers forward-orientation reading-frame configurations;
- restricts sequence modification to synonymous codon substitutions; and
- evaluates candidate CERs primarily according to nucleotide-level agreement.

The current results demonstrate computational sequence compatibility and
amino-acid sequence preservation. They do not by themselves establish successful
biological expression or function of an optimised co-encoding construct.

Future development will include evaluation across additional protein pairs,
reverse-complement and antisense configurations, additional biologically relevant
optimisation constraints, and experimental validation.

## Author

Sofie van de Bijlmer

## License
