# acp-genome-mining-workflow

A documented, web-server-based workflow for shortlisting candidate anticancer peptides (ACPs) from the biosynthetic gene clusters (BGCs) of a bacterial genome, plus a script that screens and ranks candidates. Everything here is **computational prediction only**. No peptide has been synthesised or tested in the lab.

## Status

- Related manuscript: submitted for peer review, not yet published.
- Sequences, scores and result tables will be added after the manuscript decision.
- This repository currently documents the method and the ranking logic.

## Question

Can precursor peptides encoded in the biosynthetic gene clusters of a bacterial genome be prioritised as candidate anticancer peptides, using sequence-based predictors and simple physicochemical screening?

## Workflow

| Step | Purpose | Tool |
|------|---------|------|
| 1 | Obtain a finished bacterial genome | NCBI GenBank |
| 2 | Detect BGCs and list their core biosynthetic genes | antiSMASH 8.0.4 (web server, relaxed strictness) |
| 3 | Separate short precursor peptides from large enzymes, using domain annotation and protein length | Manual classification |
| 4 | Predict antimicrobial peptide regions | CAMPR3 |
| 5 | Design and score single-residue variants for anticancer potential | AntiCP 2.0 (Design module, SVM threshold 0.0, main ACP/AMP model) |
| 6 | Calculate net charge, hydrophobic ratio and Boman index | APD3 |
| 7 | Calculate instability index, GRAVY and half-life estimate | ExPASy ProtParam |
| 8 | Screen, score and rank candidates | `scripts/rank_peptides.py` |

Steps 1 to 7 are run on public web servers. Step 8 is the only step done in code.

## Screening and ranking logic

**Screen.** A peptide is kept only if all three conditions hold:
- net charge above 0
- instability index below 40
- Boman index below 1 kcal/mol

**Composite score.**
1. Each
