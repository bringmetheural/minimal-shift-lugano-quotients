# Minimal quotients for the SHIFT/Lugano causal witness

This repository contains the manuscript, data, figures, and reproducibility code for:

**Minimal quotients for the SHIFT/Lugano causal witness**  
Dmitry Alexandrovich Lugin  
Independent researcher

Repository: https://github.com/bringmetheural/minimal-shift-lugano-quotients

DOI: https://doi.org/10.5281/zenodo.20128994

## Main result

The SHIFT measurement is an eight-outcome product-basis measurement related to the Lugano causal witness.

The paper proves that full eight-outcome resolution is not operationally necessary for the Lugano witness. Let `L` be the eight-valued SHIFT outcome and let `C = f(L)` be the four-valued Lugano output class induced by

```text
f(0)=f(1)=0,   f(+)=f(-)=1.
```

For any stochastic coarse-graining of `L` to at most `m` outcomes, the optimal Lugano decoding success obeys

```text
P_success <= min(m,4)/4.
```

Therefore:

```text
m <= 3  -> cannot violate the Lugano causal bound 3/4;
m = 4   -> the quotient Y=f(L) reaches P_success = 1.
```

The repository includes exhaustive verification over all 4140 deterministic partitions of the eight SHIFT labels.

## Repository structure

```text
paper/     LaTeX manuscript, bibliography, PDF, and figures
scripts/   Python scripts for reproducing the calculations
data/      CSV outputs used in the manuscript
proofs/    Supplementary proof notes
```

## Reproducibility

Install the Python dependencies:

```bash
pip install -r requirements.txt
```

Run the numerical reproduction script:

```bash
python scripts/reproduce_all.py --out results
```

This regenerates the core CSV tables and figures.

## Citation

Please cite the archived release using DOI: https://doi.org/10.5281/zenodo.20128994.

Please use the metadata in `CITATION.cff`.

## License

Code is released under the MIT License.
Paper text and figures are released under the Creative Commons Attribution 4.0 International License.
