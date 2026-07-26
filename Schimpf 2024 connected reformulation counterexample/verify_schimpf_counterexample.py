#!/usr/bin/env python3
"""Run the end-to-end exact connected-pole verification."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from schimpf_counterexample.connected import compute_connected_quotient
from schimpf_counterexample.paper_data import assert_matches_paper
from schimpf_counterexample.vertex import (
    assert_vertex_tables_equal,
    derive_vertex_table,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cutoff",
        type=int,
        default=18,
        help="maximum quotient length (default: 18)",
    )
    parser.add_argument(
        "--stability-cutoff",
        type=int,
        default=16,
        help="shorter cutoff used for reconstruction stability (default: 16)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="optional path for a JSON verification report",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.stability_cutoff >= args.cutoff:
        raise ValueError("stability cutoff must be smaller than cutoff")
    stability_table = derive_vertex_table(
        maximum_length=args.stability_cutoff
    )
    table = derive_vertex_table(maximum_length=args.cutoff)
    assert_vertex_tables_equal(stability_table, table)
    if args.cutoff == 18:
        assert_matches_paper(table)
    result = compute_connected_quotient(table)

    print(
        "Rational reconstruction stability "
        f"({args.stability_cutoff} vs {args.cutoff}): PASS"
    )
    print(f"C_{{1,1}}(p) = {result.C11}")
    print(f"C_{{2,1}}(p) = {result.C21}")
    print(
        "Exact target: 18*p^2/((p-1)*(p+1)^3) "
        f"{'PASS' if result.C21 == result.expected_C21 else 'FAIL'}"
    )
    print(f"Pole order at p=1: {-result.leading_order}")
    print(f"Residue at p=1: {result.residue}")
    print(
        "Higher principal-part coefficients: "
        f"(p-1)^-3={result.laurent[-3]}, "
        f"(p-1)^-2={result.laurent[-2]}"
    )
    if result.C11.is_zero() and result.residue == Fraction(9, 4):
        print("VERIFIED: the connected coefficient has residue 9/4 at p=1.")
    else:
        raise AssertionError("connected-pole verification failed")

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(result.to_dict(), indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
