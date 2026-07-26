"""Connected-quotient calculation in bidegree (2,1)."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .exact import (
    RationalFunction,
    laurent_at_one,
    leading_rational_at_one,
    multiply_polynomials,
    polynomial,
)
from .vertex import RelativePartitionPair, VertexTable


def polynomial_power(values: tuple[int, ...], exponent: int):
    result = polynomial([1])
    base = polynomial(values)
    for _ in range(exponent):
        result = multiply_polynomials(result, base)
    return result


def toric_chain_weight_table() -> dict[RelativePartitionPair, RationalFunction]:
    p_plus_one_2 = polynomial_power((1, 1), 2)
    p_plus_one_4 = polynomial_power((1, 1), 4)
    p_plus_one_6 = polynomial_power((1, 1), 6)
    p_minus_one_2 = polynomial_power((-1, 1), 2)
    return {
        ((), ()): RationalFunction([1]),
        ((1,), ()): RationalFunction([0, 1], p_plus_one_2),
        ((), (1,)): RationalFunction([0, 1], p_plus_one_2),
        ((1,), (1,)): RationalFunction(
            [0, -1, -1, -1], p_plus_one_4
        ),
        ((2,), ()): RationalFunction(
            [0, 0, 0, -1],
            multiply_polynomials(p_minus_one_2, p_plus_one_4),
        ),
        ((1, 1), ()): RationalFunction(
            [0, 0, 0, -1],
            multiply_polynomials(p_minus_one_2, p_plus_one_4),
        ),
        ((2,), (1,)): RationalFunction(
            [0, 0, -1, -1, 0, 1],
            multiply_polynomials(p_minus_one_2, p_plus_one_6),
        ),
        ((1, 1), (1,)): RationalFunction(
            [0, 0, 0, 1, 0, -1, -1],
            multiply_polynomials(p_minus_one_2, p_plus_one_6),
        ),
    }


@dataclass(frozen=True)
class LaurentTerm:
    label: str
    value: RationalFunction
    coefficients: dict[int, Fraction]

    def to_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "value": self.value.to_dict(),
            "laurent_at_p_equals_1": {
                str(power): str(coefficient)
                for power, coefficient in self.coefficients.items()
            },
        }


@dataclass(frozen=True)
class ConnectedResult:
    A: dict[RelativePartitionPair, RationalFunction]
    B: dict[RelativePartitionPair, RationalFunction]
    aggregates: dict[str, RationalFunction]
    C11: RationalFunction
    C21: RationalFunction
    expected_C21: RationalFunction
    leading_order: int
    residue: Fraction
    laurent: dict[int, Fraction]
    terms: tuple[LaurentTerm, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "aggregates": {
                key: value.to_dict()
                for key, value in self.aggregates.items()
            },
            "C11": self.C11.to_dict(),
            "C21": self.C21.to_dict(),
            "C21_expected": self.expected_C21.to_dict(),
            "leading_order_at_p_equals_1": self.leading_order,
            "residue_at_p_equals_1": str(self.residue),
            "laurent_at_p_equals_1": {
                str(power): str(coefficient)
                for power, coefficient in self.laurent.items()
            },
            "laurent_terms": [term.to_dict() for term in self.terms],
        }


def compute_connected_quotient(table: VertexTable) -> ConnectedResult:
    A = toric_chain_weight_table()
    if set(A) != set(table.by_pair):
        raise AssertionError(
            "vertex and toric-chain relative-partition pairs differ"
        )
    B = {
        pair: A[pair] * row.h
        for pair, row in table.by_pair.items()
    }

    Z10 = A[((1,), ())]
    Z01 = A[((), (1,))]
    Z11 = A[((1,), (1,))]
    Z20 = A[((2,), ())] + A[((1, 1), ())]
    F10 = B[((1,), ())]
    F01 = B[((), (1,))]
    F11 = B[((1,), (1,))]
    F20 = B[((2,), ())] + B[((1, 1), ())]
    F21 = B[((2,), (1,))] + B[((1, 1), (1,))]

    C11 = F11 - Z10 * F01 - Z01 * F10
    term_values = (
        ("F_{2,1}", F21),
        ("-Z_{0,1}F_{2,0}", -Z01 * F20),
        ("-Z_{1,0}F_{1,1}", -Z10 * F11),
        ("-Z_{2,0}F_{0,1}", -Z20 * F01),
        ("-Z_{1,1}F_{1,0}", -Z11 * F10),
        ("2Z_{1,0}Z_{0,1}F_{1,0}", 2 * Z10 * Z01 * F10),
        ("Z_{1,0}^2F_{0,1}", Z10 * Z10 * F01),
    )
    C21 = sum((value for _, value in term_values), RationalFunction([0]))
    expected_C21 = RationalFunction(
        [0, 0, 18],
        multiply_polynomials(
            polynomial([-1, 1]), polynomial_power((1, 1), 3)
        ),
    )
    if C21 != expected_C21:
        raise AssertionError(
            "C21 does not equal 18*p^2/((p-1)*(p+1)^3)"
        )

    leading_order, residue = leading_rational_at_one(C21)
    if leading_order is None:
        raise AssertionError("C21 vanished identically")
    terms = tuple(
        LaurentTerm(
            label=label,
            value=value,
            coefficients=laurent_at_one(value, -3, 0),
        )
        for label, value in term_values
    )
    return ConnectedResult(
        A=A,
        B=B,
        aggregates={
            "Z10": Z10,
            "Z01": Z01,
            "Z11": Z11,
            "Z20": Z20,
            "F10": F10,
            "F01": F01,
            "F11": F11,
            "F20": F20,
            "F21": F21,
        },
        C11=C11,
        C21=C21,
        expected_C21=expected_C21,
        leading_order=leading_order,
        residue=residue,
        laurent=laurent_at_one(C21, -3, 0),
        terms=terms,
    )
