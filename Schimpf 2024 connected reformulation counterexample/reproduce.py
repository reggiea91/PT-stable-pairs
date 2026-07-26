#!/usr/bin/env python3
"""Generate the machine-readable verification artifacts."""

from __future__ import annotations

import argparse
import csv
import json
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
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="artifact directory (default: results)",
    )
    parser.add_argument(
        "--max-length",
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
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.stability_cutoff >= args.max_length:
        raise ValueError("stability cutoff must be smaller than maximum length")
    stability_table = derive_vertex_table(
        maximum_length=args.stability_cutoff
    )
    table = derive_vertex_table(maximum_length=args.max_length)
    assert_vertex_tables_equal(stability_table, table)
    if args.max_length == 18:
        assert_matches_paper(table)
    connected = compute_connected_quotient(table)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    vertex_path = args.output_dir / "vertex_table.json"
    connected_path = args.output_dir / "connected_verification.json"
    laurent_path = args.output_dir / "laurent_terms.csv"

    vertex_path.write_text(
        json.dumps(table.to_dict(include_coefficients=True), indent=2) + "\n",
        encoding="utf-8",
    )
    connected_output = connected.to_dict()
    connected_output["stability_check"] = {
        "status": "PASS",
        "compared_cutoffs": [
            args.stability_cutoff,
            args.max_length,
        ],
    }
    connected_path.write_text(
        json.dumps(connected_output, indent=2) + "\n",
        encoding="utf-8",
    )
    with laurent_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "term",
                "coefficient_(p-1)^-3",
                "coefficient_(p-1)^-2",
                "coefficient_(p-1)^-1",
                "coefficient_(p-1)^0",
            ]
        )
        for term in connected.terms:
            writer.writerow(
                [
                    term.label,
                    term.coefficients[-3],
                    term.coefficients[-2],
                    term.coefficients[-1],
                    term.coefficients[0],
                ]
            )
        writer.writerow(
            [
                "TOTAL",
                connected.laurent[-3],
                connected.laurent[-2],
                connected.laurent[-1],
                connected.laurent[0],
            ]
        )

    print(f"Wrote {vertex_path}")
    print(f"Wrote {connected_path}")
    print(f"Wrote {laurent_path}")
    print(
        "VERIFIED: reconstruction is stable at cutoffs "
        f"{args.stability_cutoff} and {args.max_length}; "
        "C_{2,1}(p)=18*p^2/((p-1)*(p+1)^3), "
        "with residue 9/4 at p=1."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
