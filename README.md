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
1. Each property is scaled to 0-1 across the peptides that passed the screen (min-max scaling).
2. The Boman and instability indices are inverted and half-life is log-transformed, so a higher value is always better.
3. The score combines 70% membrane-activity properties (net charge, hydrophobic ratio, GRAVY, inverted Boman index) and 30% stability (inverted instability index, log half-life), with equal weight inside each group.

The AntiCP score is recorded and compared with the ranking, but it is not part of the composite score, so the two measures remain independent.

## Using the script

Requirements: Python 3, `pandas`, `numpy`, `scipy`.

```bash
pip install pandas numpy scipy
python scripts/rank_peptides.py example/example_input.csv ranked_output.csv
```

The input CSV needs these columns: `peptide_id`, `net_charge`, `hydrophobic_ratio`, `boman_index`, `gravy`, `instability_index`, `half_life_h`, `anticp_score`. The file in `example/` uses invented placeholder values to show the format. They are not real peptides or results.

## Limitations

- All outputs are predictions. Anticancer activity, toxicity and haemolysis have not been tested.
- The 70:30 weighting is a judgement call. Because scaling depends on which peptides pass the screen, small score differences should not be over-interpreted.
- Half-life is estimated from the N-terminal residue and says nothing about stability in blood serum.
- Removal of the leader sequence and chemical modification of the mature peptide are not modelled.
- Only one genome has been analysed.

## Planned additions

- Results and sequences, after the manuscript decision
- Toxicity and haemolysis prediction steps
- Application to additional genomes

## Author

Adewuyi Joel Adeyemi. Contact: joeladewuyi8@gmail.com

The ranking script was written with AI assistance. The author is reviewing it against the manuscript's results.
