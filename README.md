# Minimal quotients for the SHIFT/Lugano causal witness

This repository contains the manuscript, data, figures, and reproduction code for the paper:

**Minimal quotients for the SHIFT/Lugano causal witness**  
Dmitry Alexandrovich Lugin  
Independent researcher

Repository: <https://github.com/bringmetheural/minimal-shift-lugano-quotients>

## Main result

The paper studies coarse-grained versions of the eight-outcome SHIFT measurement used in the SHIFT/Lugano causal-witness setting.

Let `L` be the eight-valued SHIFT outcome and let `C = f(L)` be the four-valued Lugano output class induced by

```text
f(0)=f(1)=0,   f(+)=f(-)=1.
```

For any stochastic coarse-graining of `L` to at most `m` outcomes, the optimal Lugano decoding success satisfies

```text
P_success <= min(m,4)/4.
```

Consequently:

```text
m <= 3  -> cannot violate the Lugano causal bound 3/4;
m = 4   -> the quotient Y=f(L) reaches P_success = 1.
```

The repository includes an exhaustive check over all `4140` deterministic partitions of the eight SHIFT labels.

## Contents

```text
paper/      LaTeX manuscript, bibliography, PDF, and figures
data/       CSV tables used in the manuscript
scripts/    Python reproduction scripts
proofs/     Proposition proof in Markdown and LaTeX form
```

## Reproducibility

The main numerical checks and figures can be regenerated with:

```bash
python scripts/reproduce_all.py --out results
```

The script writes regenerated tables and figures to `results/`.

Python dependencies are listed in `requirements.txt`.

## Citation

Citation metadata is provided in `CITATION.cff`.

## License

- Code: MIT License, see `LICENSE`.
- Paper text and figures: Creative Commons Attribution 4.0 International, see `LICENSE-CC-BY-4.0.txt`.
