"""Derivation of the eight vertex rows used in bidegree (2,1)."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .boxes import (
    Partition,
    enumerate_relative_partition_pair_series,
    reconstruct_rational_series,
)
from .exact import RationalFunction

RelativePartitionPair = tuple[Partition, Partition]

VERTEX_RELATIVE_PARTITION_PAIRS: tuple[RelativePartitionPair, ...] = (
    ((), ()),
    ((1,), ()),
    ((), (1,)),
    ((1,), (1,)),
    ((2,), ()),
    ((1, 1), ()),
    ((2,), (1,)),
    ((1, 1), (1,)),
)


def partition_label(partition: Partition) -> str:
    return r"\varnothing" if not partition else (
        "(" + ",".join(str(part) for part in partition) + ")"
    )


def relative_partition_pair_label(pair: RelativePartitionPair) -> str:
    lambda_partition, mu_partition = pair
    return (
        f"{partition_label(lambda_partition)}"
        f"|{partition_label(mu_partition)}"
    )


@dataclass(frozen=True)
class VertexRow:
    pair: RelativePartitionPair
    configuration_count: int
    P_coefficients: tuple[Fraction, ...]
    N_coefficients: tuple[Fraction, ...]
    P: RationalFunction
    N: RationalFunction
    h: RationalFunction

    @property
    def label(self) -> str:
        return relative_partition_pair_label(self.pair)

    def to_dict(self, include_coefficients: bool = False) -> dict[str, object]:
        lambda_partition, mu_partition = self.pair
        output: dict[str, object] = {
            "label": self.label,
            "lambda": list(lambda_partition),
            "mu": list(mu_partition),
            "configuration_count": self.configuration_count,
            "P": self.P.to_dict(),
            "N": self.N.to_dict(),
            "h": self.h.to_dict(),
        }
        if include_coefficients:
            output["enumerated_coefficients"] = {
                "P": [str(value) for value in self.P_coefficients],
                "N": [str(value) for value in self.N_coefficients],
            }
        return output


@dataclass(frozen=True)
class VertexTable:
    maximum_length: int
    degree: int
    weights: tuple[Fraction, Fraction, Fraction]
    rows: tuple[VertexRow, ...]

    @property
    def by_pair(self) -> dict[RelativePartitionPair, VertexRow]:
        return {row.pair: row for row in self.rows}

    def to_dict(self, include_coefficients: bool = False) -> dict[str, object]:
        return {
            "maximum_length": self.maximum_length,
            "degree": self.degree,
            "weights": [str(value) for value in self.weights],
            "rows": [
                row.to_dict(include_coefficients=include_coefficients)
                for row in self.rows
            ],
        }


def assert_vertex_tables_equal(left: VertexTable, right: VertexTable) -> None:
    if left.degree != right.degree:
        raise AssertionError("vertex tables use different descendent degrees")
    if left.weights != right.weights:
        raise AssertionError("vertex tables use different torus weights")
    if set(left.by_pair) != set(right.by_pair):
        raise AssertionError(
            "vertex tables have different relative-partition pairs"
        )
    for pair, left_row in left.by_pair.items():
        right_row = right.by_pair[pair]
        for name in ("P", "N", "h"):
            if getattr(left_row, name) != getattr(right_row, name):
                raise AssertionError(
                    f"{name} is not stable for {left_row.label}"
                )


def derive_vertex_table(
    maximum_length: int = 18,
    degree: int = 3,
    weights: tuple[Fraction, Fraction, Fraction] = (
        Fraction(1),
        Fraction(2),
        Fraction(-3),
    ),
) -> VertexTable:
    rows: list[VertexRow] = []
    for pair in VERTEX_RELATIVE_PARTITION_PAIRS:
        lambda_partition, mu_partition = pair
        series = enumerate_relative_partition_pair_series(
            lambda_partition,
            mu_partition,
            maximum_length,
            degree,
            weights,
        )
        P_coefficients = series["P_coefficients"]
        N_coefficients = series["N_coefficients"]
        if not isinstance(P_coefficients, tuple) or not isinstance(
            N_coefficients, tuple
        ):
            raise TypeError("enumerated coefficients must be tuples")
        P = reconstruct_rational_series(
            P_coefficients,
            maximum_denominator_degree=14,
            maximum_numerator_degree=10,
            guard=5,
        )
        N = reconstruct_rational_series(
            N_coefficients,
            maximum_denominator_degree=18,
            maximum_numerator_degree=12,
            guard=5,
        )
        rows.append(
            VertexRow(
                pair=pair,
                configuration_count=int(series["configuration_count"]),
                P_coefficients=P_coefficients,
                N_coefficients=N_coefficients,
                P=P,
                N=N,
                h=N / P,
            )
        )
    return VertexTable(
        maximum_length=maximum_length,
        degree=degree,
        weights=weights,
        rows=tuple(rows),
    )
