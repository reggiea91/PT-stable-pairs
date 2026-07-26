# Exact verification of a connected stable pairs pole counterexample

This repository contains independent Python and JavaScript implementations of
the finite calculation in:

> Reginald Anderson, *A toric Calabi--Yau counterexample to Schimpf's connected stable pairs pole conjecture*.

For the curve class

\[
\beta=2[C_1]+[C_2],\qquad \operatorname{div}(\beta)=1,
\]

the calculation gives the connected coefficient

\[
C_{2,1}(p)=\frac{18p^2}{(p-1)(p+1)^3}.
\]

It therefore has a simple pole at \(p=1\) with residue \(9/4\).

The mathematical proof is self-contained in the paper. This repository
independently reproduces the finite two-leg vertex enumeration and the exact
rational-function calculations used there.

## Choose one reproduction route

The Python and JavaScript implementations are alternative, complete
reproduction routes. A reader may choose either one; it is not necessary to
run both. Both routes use the same length cutoffs \(16\) and \(18\), reconstruct
the same eight vertex rows, form the same connected quotient, and verify the
same Laurent expansion at \(p=1\). Running the second route is an optional
independent cross-check.

The paper does not need to be compiled, and no TeX source is required.

## Option A: Python

Requirement: Python 3.10 or newer. No third-party Python packages are needed.

From the repository root, run:

```console
python -m unittest discover -s tests -v
python verify_schimpf_counterexample.py
python reproduce.py
```

The test suite reports five passing tests. The verifier concludes with:

```text
Rational reconstruction stability (16 vs 18): PASS
C_{1,1}(p) = 0
C_{2,1}(p) = (-18*p^2)/(1 + 2*p - 2*p^3 - p^4)
Exact target: 18*p^2/((p-1)*(p+1)^3) PASS
Pole order at p=1: 1
Residue at p=1: 9/4
Higher principal-part coefficients: (p-1)^-3=0, (p-1)^-2=0
VERIFIED: the connected coefficient has residue 9/4 at p=1.
```

To print the eight reconstructed vertex rows separately, run:

```console
python pt_two_leg_box_enum.py
```

## Option B: JavaScript

Requirements:

- Node.js 20 or newer
- npm

From the repository root, run:

```console
npm ci
npm run check
```

The test suite reports five passing tests. The command then generates the exact
certificate under `results/` and concludes with:

```text
Finite two-leg PT vertex enumeration: PASS
Rational reconstruction stability (16 vs 18): PASS
C11(p) = 0
C21(p) = (-18*p^2)/(1 + 2*p - 2*p^3 - p^4)
Exact identity C21(p) = 18*p^2/((p-1)*(p+1)^3): PASS
Residue at p=1: 9/4
Forbidden simple pole at p=1: VERIFIED
```

## Mathematical pipeline

Whichever route is chosen performs the following complete calculation:

1. Enumerate quotient-box configurations through lengths \(16\) and \(18\) for
   the eight pairs of relative partitions needed in bidegree \((2,1)\).
2. Form the no-descendent and \([z^3]\)-descendent coefficient series.
3. Reconstruct the exact rational functions
   \(P_{\lambda\mid\mu}\) and \(N_{\lambda\mid\mu}\).
4. Verify that every reconstructed rational function agrees at cutoffs \(16\)
   and \(18\).
5. Compare all eight reconstructed rows with the formulas in the paper.
6. Compute
   \(h_{\lambda\mid\mu}=N_{\lambda\mid\mu}/P_{\lambda\mid\mu}\).
7. Combine these ratios with the explicitly recorded toric-chain contributions
   \(A_{\lambda\mid\mu}\), and form \(Z_{a,b}\) and \(F_{a,b}\).
8. Form the connected quotient coefficient
   \[
   \begin{aligned}
   C_{2,1}={}&F_{2,1}-Z_{0,1}F_{2,0}-Z_{1,0}F_{1,1}
   -Z_{2,0}F_{0,1}-Z_{1,1}F_{1,0}\\
   &+2Z_{1,0}Z_{0,1}F_{1,0}+Z_{1,0}^2F_{0,1}.
   \end{aligned}
   \]
9. Simplify the result exactly and extract its Laurent coefficients at
   \(p=1\).

The implementations do not call one another and do not transfer computed
rational functions between languages.

## Generated results

Option B writes:

- `results/certificate.json`
- `results/laurent-ledger.csv`
- `results/SHA256SUMS`

Option A writes:

- `results/vertex_table.json`
- `results/connected_verification.json`
- `results/laurent_terms.csv`

The JSON files contain the reconstructed vertex data, bidegree aggregates,
connected quotient, Laurent coefficients, and exact conclusion. The CSV files
give the term-by-term Laurent calculation.

## Repository structure

- `src/` and `scripts/`: JavaScript implementation
- `test/`: JavaScript tests
- `schimpf_counterexample/`: Python implementation
- `tests/`: Python tests
- `results/`: exact generated certificates
- `.github/workflows/`: continuous verification
- `docs/REPRODUCIBILITY.md`: detailed reproduction information

## Citation

Citation metadata is provided in `CITATION.cff`. The archived release
corresponding to the paper should be cited by its version-specific Zenodo DOI.
