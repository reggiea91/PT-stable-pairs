#!/usr/bin/env python3
"""Reconstruct the printed two-leg vertex table from finite configurations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from schimpf_counterexample.paper_data import assert_matches_paper
from schimpf_counterexample.vertex import derive_vertex_table


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cutoff",
        type=int,
        default=18,
        help="maximum quotient length (default: 18)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="optional path for a JSON copy of the derived table",
    )
    parser.add_argument(
        "--include-coefficients",
        action="store_true",
        help="include all enumerated coefficients in JSON output",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    table = derive_vertex_table(maximum_length=args.cutoff)
    if args.cutoff == 18:
        assert_matches_paper(table)
    print(
        "Finite two-leg vertex reconstruction "
        f"(cutoff={table.maximum_length}, degree={table.degree}, "
        "weights=(1,2,-3))"
    )
    for row in table.rows:
        print(f"\n{row.label}; configurations={row.configuration_count}")
        print(f"  P = {row.P}")
        print(f"  N = {row.N}")
        print(f"  h = N/P = {row.h}")
    if args.cutoff == 18:
        print("\nPASS: every row agrees exactly with the accompanying paper.")

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(
                table.to_dict(
                    include_coefficients=args.include_coefficients
                ),
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
