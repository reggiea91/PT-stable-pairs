"""Verification tests for the exact counterexample computation."""

from __future__ import annotations

import unittest
from fractions import Fraction

from schimpf_counterexample.connected import compute_connected_quotient
from schimpf_counterexample.exact import (
    RationalFunction,
    laurent_at_one,
    leading_rational_at_one,
)
from schimpf_counterexample.paper_data import (
    PAPER_VERTEX_DATA,
    assert_matches_paper,
)
from schimpf_counterexample.vertex import (
    assert_vertex_tables_equal,
    derive_vertex_table,
)


class ExactArithmeticTests(unittest.TestCase):
    def test_reduction_and_laurent_expansion(self) -> None:
        value = RationalFunction(
            [0, 0, 18], [-1, -2, 0, 2, 1]
        )
        order, coefficient = leading_rational_at_one(value)
        self.assertEqual(order, -1)
        self.assertEqual(coefficient, Fraction(9, 4))
        self.assertEqual(
            laurent_at_one(value, -1, 0),
            {-1: Fraction(9, 4), 0: Fraction(9, 8)},
        )


class CounterexampleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.table = derive_vertex_table(maximum_length=18)
        cls.result = compute_connected_quotient(cls.table)

    def test_all_vertex_rows_match_the_paper(self) -> None:
        assert_matches_paper(self.table)
        self.assertEqual(len(self.table.rows), 8)
        self.assertEqual(set(self.table.by_pair), set(PAPER_VERTEX_DATA))

    def test_vertex_reconstruction_stabilizes(self) -> None:
        shorter = derive_vertex_table(maximum_length=16)
        assert_vertex_tables_equal(shorter, self.table)

    def test_connected_identity_and_residue(self) -> None:
        self.assertTrue(self.result.C11.is_zero())
        self.assertEqual(self.result.C21, self.result.expected_C21)
        self.assertEqual(self.result.leading_order, -1)
        self.assertEqual(self.result.residue, Fraction(9, 4))
        self.assertEqual(self.result.laurent[-3], 0)
        self.assertEqual(self.result.laurent[-2], 0)
        self.assertEqual(self.result.laurent[-1], Fraction(9, 4))
        self.assertEqual(self.result.laurent[0], Fraction(9, 8))

    def test_printed_term_by_term_laurent_table(self) -> None:
        expected = {
            "F_{2,1}": (
                Fraction(-3, 16),
                Fraction(-9, 32),
                Fraction(69, 32),
            ),
            "-Z_{0,1}F_{2,0}": (
                Fraction(3, 16),
                Fraction(9, 32),
                Fraction(3, 64),
            ),
            "-Z_{1,0}F_{1,1}": (Fraction(0),) * 3,
            "-Z_{2,0}F_{0,1}": (
                Fraction(0),
                Fraction(0),
                Fraction(3, 64),
            ),
            "-Z_{1,1}F_{1,0}": (Fraction(0),) * 3,
            "2Z_{1,0}Z_{0,1}F_{1,0}": (Fraction(0),) * 3,
            "Z_{1,0}^2F_{0,1}": (Fraction(0),) * 3,
        }
        for term in self.result.terms:
            self.assertEqual(
                (
                    term.coefficients[-3],
                    term.coefficients[-2],
                    term.coefficients[-1],
                ),
                expected[term.label],
            )

if __name__ == "__main__":
    unittest.main()
