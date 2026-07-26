"""Exact arithmetic for rational functions in one variable.

Polynomials are tuples of ``fractions.Fraction`` coefficients in ascending
order.  Thus ``(a_0, a_1, ...)`` represents ``a_0 + a_1*p + ...``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb
from typing import Iterable

Polynomial = tuple[Fraction, ...]


def polynomial(values: Iterable[int | Fraction]) -> Polynomial:
    return trim(tuple(Fraction(value) for value in values))


def trim(values: Polynomial) -> Polynomial:
    if not values:
        return (Fraction(0),)
    length = len(values)
    while length > 1 and values[length - 1] == 0:
        length -= 1
    return values[:length]


ZERO_POLYNOMIAL = polynomial([0])
ONE_POLYNOMIAL = polynomial([1])


def is_zero(values: Polynomial) -> bool:
    return trim(values) == ZERO_POLYNOMIAL


def add_polynomials(left: Polynomial, right: Polynomial) -> Polynomial:
    length = max(len(left), len(right))
    return trim(
        tuple(
            (left[index] if index < len(left) else Fraction(0))
            + (right[index] if index < len(right) else Fraction(0))
            for index in range(length)
        )
    )


def negate_polynomial(values: Polynomial) -> Polynomial:
    return tuple(-coefficient for coefficient in values)


def subtract_polynomials(left: Polynomial, right: Polynomial) -> Polynomial:
    return add_polynomials(left, negate_polynomial(right))


def scale_polynomial(
    values: Polynomial, scalar: int | Fraction
) -> Polynomial:
    factor = Fraction(scalar)
    return trim(tuple(coefficient * factor for coefficient in values))


def multiply_polynomials(left: Polynomial, right: Polynomial) -> Polynomial:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            result[left_degree + right_degree] += (
                left_coefficient * right_coefficient
            )
    return trim(tuple(result))


def divide_polynomials(
    dividend: Polynomial, divisor: Polynomial
) -> tuple[Polynomial, Polynomial]:
    denominator = trim(divisor)
    if is_zero(denominator):
        raise ZeroDivisionError("polynomial division by zero")
    remainder = trim(dividend)
    quotient = [Fraction(0)] * max(
        1, len(remainder) - len(denominator) + 1
    )
    while not is_zero(remainder) and len(remainder) >= len(denominator):
        degree = len(remainder) - len(denominator)
        factor = remainder[-1] / denominator[-1]
        quotient[degree] += factor
        shifted = (Fraction(0),) * degree + scale_polynomial(
            denominator, factor
        )
        remainder = subtract_polynomials(remainder, shifted)
    return trim(tuple(quotient)), trim(remainder)


def polynomial_gcd(left: Polynomial, right: Polynomial) -> Polynomial:
    a = trim(left)
    b = trim(right)
    while not is_zero(b):
        _, remainder = divide_polynomials(a, b)
        a, b = b, remainder
    if is_zero(a):
        return ONE_POLYNOMIAL
    return scale_polynomial(a, Fraction(1, 1) / a[-1])


def evaluate_polynomial(
    values: Polynomial, point: int | Fraction
) -> Fraction:
    x = Fraction(point)
    result = Fraction(0)
    for coefficient in reversed(values):
        result = result * x + coefficient
    return result


@dataclass(frozen=True, init=False)
class RationalFunction:
    """A reduced rational function over the rational numbers."""

    numerator: Polynomial
    denominator: Polynomial

    def __init__(
        self,
        numerator: Iterable[int | Fraction],
        denominator: Iterable[int | Fraction] = (1,),
    ) -> None:
        top = polynomial(numerator)
        bottom = polynomial(denominator)
        if is_zero(bottom):
            raise ZeroDivisionError("zero rational-function denominator")
        if is_zero(top):
            object.__setattr__(self, "numerator", ZERO_POLYNOMIAL)
            object.__setattr__(self, "denominator", ONE_POLYNOMIAL)
            return

        divisor = polynomial_gcd(top, bottom)
        top, top_remainder = divide_polynomials(top, divisor)
        bottom, bottom_remainder = divide_polynomials(bottom, divisor)
        if not is_zero(top_remainder) or not is_zero(bottom_remainder):
            raise ArithmeticError("polynomial gcd did not divide exactly")

        normalization = bottom[0] if bottom[0] != 0 else bottom[-1]
        top = scale_polynomial(top, Fraction(1, 1) / normalization)
        bottom = scale_polynomial(bottom, Fraction(1, 1) / normalization)
        object.__setattr__(self, "numerator", top)
        object.__setattr__(self, "denominator", bottom)

    @staticmethod
    def coerce(value: RationalFunction | int | Fraction) -> RationalFunction:
        if isinstance(value, RationalFunction):
            return value
        return RationalFunction([Fraction(value)])

    def __add__(
        self, other: RationalFunction | int | Fraction
    ) -> RationalFunction:
        right = RationalFunction.coerce(other)
        return RationalFunction(
            add_polynomials(
                multiply_polynomials(self.numerator, right.denominator),
                multiply_polynomials(right.numerator, self.denominator),
            ),
            multiply_polynomials(self.denominator, right.denominator),
        )

    def __radd__(
        self, other: RationalFunction | int | Fraction
    ) -> RationalFunction:
        return self + other

    def __neg__(self) -> RationalFunction:
        return RationalFunction(
            negate_polynomial(self.numerator), self.denominator
        )

    def __sub__(
        self, other: RationalFunction | int | Fraction
    ) -> RationalFunction:
        return self + (-RationalFunction.coerce(other))

    def __rsub__(
        self, other: RationalFunction | int | Fraction
    ) -> RationalFunction:
        return RationalFunction.coerce(other) - self

    def __mul__(
        self, other: RationalFunction | int | Fraction
    ) -> RationalFunction:
        right = RationalFunction.coerce(other)
        return RationalFunction(
            multiply_polynomials(self.numerator, right.numerator),
            multiply_polynomials(self.denominator, right.denominator),
        )

    def __rmul__(
        self, other: RationalFunction | int | Fraction
    ) -> RationalFunction:
        return self * other

    def __truediv__(
        self, other: RationalFunction | int | Fraction
    ) -> RationalFunction:
        right = RationalFunction.coerce(other)
        if is_zero(right.numerator):
            raise ZeroDivisionError("rational-function division by zero")
        return RationalFunction(
            multiply_polynomials(self.numerator, right.denominator),
            multiply_polynomials(self.denominator, right.numerator),
        )

    def __rtruediv__(
        self, other: RationalFunction | int | Fraction
    ) -> RationalFunction:
        return RationalFunction.coerce(other) / self

    def is_zero(self) -> bool:
        return is_zero(self.numerator)

    def to_dict(self) -> dict[str, object]:
        return {
            "numerator": [str(value) for value in self.numerator],
            "denominator": [str(value) for value in self.denominator],
            "expression": str(self),
        }

    def __str__(self) -> str:
        top = polynomial_to_string(self.numerator)
        bottom = polynomial_to_string(self.denominator)
        return top if bottom == "1" else f"({top})/({bottom})"


ZERO = RationalFunction([0])
ONE = RationalFunction([1])


def shifted_at_one(values: Polynomial) -> Polynomial:
    result = [Fraction(0)] * len(values)
    for degree, coefficient in enumerate(values):
        for index in range(degree + 1):
            result[index] += coefficient * comb(degree, index)
    return trim(tuple(result))


def leading_polynomial_at_one(
    values: Polynomial,
) -> tuple[int | None, Fraction]:
    shifted = shifted_at_one(values)
    for order, coefficient in enumerate(shifted):
        if coefficient != 0:
            return order, coefficient
    return None, Fraction(0)


def leading_rational_at_one(
    value: RationalFunction,
) -> tuple[int | None, Fraction]:
    top_order, top_coefficient = leading_polynomial_at_one(value.numerator)
    bottom_order, bottom_coefficient = leading_polynomial_at_one(
        value.denominator
    )
    if top_order is None:
        return None, Fraction(0)
    if bottom_order is None:
        raise ZeroDivisionError("identically zero denominator")
    return top_order - bottom_order, top_coefficient / bottom_coefficient


def laurent_at_one(
    value: RationalFunction, minimum_power: int, maximum_power: int
) -> dict[int, Fraction]:
    top_order, _ = leading_polynomial_at_one(value.numerator)
    bottom_order, _ = leading_polynomial_at_one(value.denominator)
    if top_order is None:
        return {
            power: Fraction(0)
            for power in range(minimum_power, maximum_power + 1)
        }
    if bottom_order is None:
        raise ZeroDivisionError("identically zero denominator")

    numerator = shifted_at_one(value.numerator)[top_order:]
    denominator = shifted_at_one(value.denominator)[bottom_order:]
    start = top_order - bottom_order
    needed = max(0, maximum_power - start)
    coefficients: list[Fraction] = []
    for index in range(needed + 1):
        right = (
            numerator[index] if index < len(numerator) else Fraction(0)
        )
        for offset in range(1, index + 1):
            denominator_coefficient = (
                denominator[offset]
                if offset < len(denominator)
                else Fraction(0)
            )
            right -= denominator_coefficient * coefficients[index - offset]
        coefficients.append(right / denominator[0])

    output: dict[int, Fraction] = {}
    for power in range(minimum_power, maximum_power + 1):
        index = power - start
        output[power] = (
            coefficients[index]
            if 0 <= index < len(coefficients)
            else Fraction(0)
        )
    return output


def polynomial_to_string(values: Polynomial, variable: str = "p") -> str:
    terms: list[tuple[str, str]] = []
    for degree, coefficient in enumerate(values):
        if coefficient == 0:
            continue
        sign = "-" if coefficient < 0 else "+"
        magnitude = abs(coefficient)
        if degree == 0:
            body = str(magnitude)
        else:
            variable_part = variable if degree == 1 else f"{variable}^{degree}"
            body = (
                variable_part
                if magnitude == 1
                else f"{magnitude}*{variable_part}"
            )
        terms.append((sign, body))
    if not terms:
        return "0"
    first_sign, first_body = terms[0]
    output = f"-{first_body}" if first_sign == "-" else first_body
    for sign, body in terms[1:]:
        output += f" {sign} {body}"
    return output
