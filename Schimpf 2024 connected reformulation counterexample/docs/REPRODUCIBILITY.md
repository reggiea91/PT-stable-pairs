# Reproducibility

## Scope

The paper contains the complete mathematical proof. The software in this
repository independently reconstructs the finite two-leg vertex table and
verifies the connected quotient coefficient

\[
C_{2,1}(p)=\frac{18p^2}{(p-1)(p+1)^3}.
\]

All arithmetic is exact.

## Choosing an implementation

Python and JavaScript provide alternative, complete reproduction routes. A
reader may use either implementation; running both is not required. Each route
uses cutoffs 16 and 18, reconstructs and checks all eight vertex rows, forms the
connected quotient, and verifies the residue \(9/4\). Running both routes gives
an optional independent cross-check.

## Requirements

The JavaScript implementation requires Node.js 20 or newer and npm. The Python
implementation requires Python 3.10 or newer. Neither implementation has
third-party runtime dependencies.

## Option A: Python

From the repository root:

```console
python -m unittest discover -s tests -v
python verify_schimpf_counterexample.py
python reproduce.py
```

The five tests compare every reconstructed vertex row with the table printed in
the paper, verify reconstruction stability at cutoffs 16 and 18, check the
connected identity, and check the term-by-term Laurent coefficients. The
verifier prints the exact coefficient and residue. The reproduction script
regenerates the Python JSON and CSV files under `results/`.

The vertex table may also be printed directly:

```console
python pt_two_leg_box_enum.py
```

## Option B: JavaScript

From the repository root:

```console
npm ci
npm run check
```

The first command installs the locked development environment. The second
command:

1. enumerates the required finite box configurations;
2. reconstructs the exact vertex table;
3. verifies stability of rational reconstruction at cutoffs 16 and 18;
4. forms the bidegree aggregates and connected quotient;
5. verifies the exact rational identity and Laurent residue;
6. regenerates the JavaScript certificate files; and
7. compares the regenerated certificate with the committed snapshot.

GitHub Actions runs the same verification on every push and pull request.

## Independent implementations

The two implementations use separate exact-arithmetic, enumeration, vertex,
and connected-quotient modules:

```text
JavaScript                         Python
----------                         ------
src/exact.js                       schimpf_counterexample/exact.py
src/box-enumerator.js              schimpf_counterexample/boxes.py
src/vertex-table.js                schimpf_counterexample/vertex.py
src/connected-quotient.js          schimpf_counterexample/connected.py
src/pipeline.js                    reproduce.py
```

Each implementation performs the same mathematical stages. The toric-chain
contributions \(A_{\lambda\mid\mu}\) are recorded explicitly in each connected
quotient module, while the \(P_{\lambda\mid\mu}\) and
\(N_{\lambda\mid\mu}\) functions are reconstructed from the enumerated
configurations. Neither implementation calls the other.

## Certificate contents

The generated files record:

- the finite-enumeration parameters;
- the reconstructed \(P_{\lambda\mid\mu}\),
  \(N_{\lambda\mid\mu}\), and \(h_{\lambda\mid\mu}\);
- the toric-chain weights and bidegree aggregates;
- the seven summands in the connected quotient;
- their Laurent coefficients at \(p=1\);
- the exact identity for \(C_{2,1}(p)\); and
- the residue \(9/4\).

`results/SHA256SUMS` gives checksums for the JavaScript-generated certificate
files.
