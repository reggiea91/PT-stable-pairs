"""Exact vertex data printed in the accompanying paper.

These constants are used only as independent verification targets.  The main
verification first derives the same functions from finite box configurations.
"""

from __future__ import annotations

from .exact import RationalFunction
from .vertex import RelativePartitionPair


def rf(
    numerator: tuple[int, ...],
    denominator: tuple[int, ...] = (1,),
) -> RationalFunction:
    return RationalFunction(numerator, denominator)


PAPER_VERTEX_DATA: dict[
    RelativePartitionPair, tuple[int, RationalFunction, RationalFunction]
] = {
    ((), ()): (1, rf((1,)), rf((0,))),
    ((1,), ()): (
        19,
        rf((1,), (1, 1)),
        rf((-3, 3), (1, 2, 1)),
    ),
    ((), (1,)): (
        19,
        rf((1,), (1, 1)),
        rf((-3, 3), (1, 2, 1)),
    ),
    ((1,), (1,)): (
        172,
        rf((1, 1, 1), (1, 2, 1)),
        rf((0, 6, -6), (1, 3, 3, 1)),
    ),
    ((2,), ()): (
        100,
        rf((1,), (1, 1, -1, -1)),
        rf((6, 6, -24), (1, 1, -2, -2, 1, 1)),
    ),
    ((1, 1), ()): (
        100,
        rf((1,), (1, 1, -1, -1)),
        rf((-24, 6, 6), (1, 1, -2, -2, 1, 1)),
    ),
    ((2,), (1,)): (
        544,
        rf((1, 1, 0, -1), (1, 2, 0, -2, -1)),
        rf(
            (15, 21, -27, -30, -12, 21),
            (1, 2, -1, -4, -1, 2, 1),
        ),
    ),
    ((1, 1), (1,)): (
        634,
        rf((1, 0, -1, -1), (1, 2, 0, -2, -1)),
        rf(
            (-21, 12, 30, 27, -21, -15),
            (1, 2, -1, -4, -1, 2, 1),
        ),
    ),
}


def assert_matches_paper(table: object) -> None:
    from .vertex import VertexTable

    if not isinstance(table, VertexTable):
        raise TypeError("expected a VertexTable")
    if set(table.by_pair) != set(PAPER_VERTEX_DATA):
        raise AssertionError("derived and printed relative-partition pairs differ")
    for pair, row in table.by_pair.items():
        count, expected_P, expected_N = PAPER_VERTEX_DATA[pair]
        if row.configuration_count != count:
            raise AssertionError(
                f"configuration count differs for {row.label}: "
                f"{row.configuration_count} != {count}"
            )
        if row.P != expected_P:
            raise AssertionError(
                f"P differs from the paper for {row.label}"
            )
        if row.N != expected_N:
            raise AssertionError(
                f"N differs from the paper for {row.label}"
            )
        if row.h != expected_N / expected_P:
            raise AssertionError(
                f"h=N/P differs from the paper for {row.label}"
            )
