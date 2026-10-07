# Crash insurance and the G10 carry premium

Amanjeet Singh. Paper, version of 7 October 2026.

- [Paper (PDF)](Crash_Insurance_and_the_G10_Carry_Premium.pdf)
- [LaTeX source](main.tex)
- [Extended abstract](summary.md)

The paper asks whether the price of crash insurance in G10 currency option smiles accounts for the carry premium and its change between the zero-rate regime (May 2013 to December 2021) and the hiking regime (from January 2022). Its estimands, samples and inference were fixed in the [research design](../research_design.md) before estimation; later changes are recorded in the [research log](../research_log.md), and post hoc analyses are labelled as such in the text. The proofs of the theory results it uses are in the [theory notes](../theory/README.md).

## Building the PDF

From this folder, with a TeX distribution that provides pdfLaTeX:

```bash
pdflatex -jobname=Crash_Insurance_and_the_G10_Carry_Premium main.tex
pdflatex -jobname=Crash_Insurance_and_the_G10_Carry_Premium main.tex
```

## Reproduction and data

The results that the [published series](../data/public/fx_carry_portfolio_series/README.md) support can be recomputed without a data licence with the command in the [reproduction instructions](../docs/reproducing_paper.md), which list them: E1 under the market reading with 10-delta hedges, E2 with φ as predictor, the secondary moment predictors, E3 and its split by rate regime, E5, the hedge-cost ratio and its confidence sets, the long at-the-money sample, the post hoc regime comparison, the comparison with published currency portfolios and part of the robustness grid. Recomputing the results that need the quotes, among them E4, the 25-delta skew prices of E1, every result of the oriented 10-delta risk reversal, the other robustness variants, the three-month tenor and the per-currency analyses, needs the reader's own LSEG Workspace licence; [licensed reproduction](../replication/licensed_reproduction.md) lists the steps, their inputs and outputs.

## Licence

The paper and its LaTeX source are released under [CC BY 4.0](../LICENSES/CC-BY-4.0.txt) within the [licensing scope](../LICENSING.md). They contain no LSEG observation: the reported numbers are pooled estimates, test statistics, intervals and coverage facts.
