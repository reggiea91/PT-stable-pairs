"""Finite two-leg box configurations and rational-series reconstruction."""

from __future__ import annotations

from fractions import Fraction
from math import factorial
from typing import Iterable

from .exact import RationalFunction

Partition = tuple[int, ...]
Box = tuple[int, int, int, int]
Weight = tuple[int, int, int]
Character = dict[Weight, int]

CUBE_TERMS = (
    (0, 0, 0, 1),
    (1, 0, 0, -1),
    (0, 1, 0, -1),
    (0, 0, 1, -1),
    (1, 1, 0, 1),
    (1, 0, 1, 1),
    (0, 1, 1, 1),
    (1, 1, 1, -1),
)
YZ_TERMS = (
    (0, 0, 0, 1),
    (0, 1, 0, -1),
    (0, 0, 1, -1),
    (0, 1, 1, 1),
)
XZ_TERMS = (
    (0, 0, 0, 1),
    (1, 0, 0, -1),
    (0, 0, 1, -1),
    (1, 0, 1, 1),
)


def partition_boxes(partition: Partition) -> tuple[tuple[int, int], ...]:
    return tuple(
        (column, row)
        for row, row_length in enumerate(partition)
        for column in range(row_length)
    )


def memberships(
    weight: Weight, lambda_partition: Partition, mu_partition: Partition
) -> tuple[int, ...]:
    i, j, k = weight
    at_weight: list[int] = []
    if (j, k) in partition_boxes(lambda_partition):
        at_weight.append(1)
    if (i, k) in partition_boxes(mu_partition):
        at_weight.append(2)
    return tuple(at_weight)


def boxes_at_weight(
    weight: Weight, lambda_partition: Partition, mu_partition: Partition
) -> tuple[Box, ...]:
    at_weight = memberships(weight, lambda_partition, mu_partition)
    if any(coordinate < 0 for coordinate in weight):
        return (
            ((at_weight[0], *weight),)
            if len(at_weight) == 1
            else ()
        )
    return ((12, *weight),) if len(at_weight) == 2 else ()


def quotient_boxes(
    lambda_partition: Partition,
    mu_partition: Partition,
    maximum_length: int,
) -> tuple[Box, ...]:
    boxes: dict[Box, Box] = {}
    lambda_boxes = partition_boxes(lambda_partition)
    mu_boxes = partition_boxes(mu_partition)
    for j, k in lambda_boxes:
        for i in range(-maximum_length, 0):
            box = (1, i, j, k)
            boxes[box] = box
    for i, k in mu_boxes:
        for j in range(-maximum_length, 0):
            box = (2, i, j, k)
            boxes[box] = box
    for j, lambda_row in lambda_boxes:
        for i, mu_row in mu_boxes:
            if lambda_row == mu_row:
                box = (12, i, j, lambda_row)
                boxes[box] = box
    return tuple(boxes.values())


def successors(
    box: Box, lambda_partition: Partition, mu_partition: Partition
) -> tuple[Box, ...]:
    _, i, j, k = box
    output: list[Box] = []
    for di, dj, dk in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        output.extend(
            boxes_at_weight(
                (i + di, j + dj, k + dk),
                lambda_partition,
                mu_partition,
            )
        )
    return tuple(output)


def enumerate_configurations(
    lambda_partition: Partition,
    mu_partition: Partition,
    maximum_length: int,
) -> tuple[frozenset[Box], ...]:
    """Enumerate upward-closed finite quotient-box configurations."""

    boxes = quotient_boxes(
        lambda_partition, mu_partition, maximum_length
    )
    available = set(boxes)
    successor_map = {
        box: tuple(
            successor
            for successor in successors(
                box, lambda_partition, mu_partition
            )
            if successor in available
        )
        for box in boxes
    }
    ordered = sorted(
        boxes,
        key=lambda box: (box[1] + box[2] + box[3], box),
        reverse=True,
    )
    chosen: set[Box] = set()
    configurations: list[frozenset[Box]] = []

    def visit(index: int) -> None:
        if len(chosen) > maximum_length:
            return
        if index == len(ordered):
            configurations.append(frozenset(chosen))
            return
        box = ordered[index]
        visit(index + 1)
        if all(successor in chosen for successor in successor_map[box]):
            chosen.add(box)
            visit(index + 1)
            chosen.remove(box)

    visit(0)
    return tuple(configurations)


def add_character(
    character: Character, weight: Weight, coefficient: int
) -> None:
    next_coefficient = character.get(weight, 0) + coefficient
    if next_coefficient == 0:
        character.pop(weight, None)
    else:
        character[weight] = next_coefficient


def add_terms(
    character: Character,
    base: Weight,
    terms: Iterable[tuple[int, int, int, int]],
    coefficient: int,
) -> None:
    for di, dj, dk, sign in terms:
        add_character(
            character,
            (base[0] + di, base[1] + dj, base[2] + dk),
            coefficient * sign,
        )


def cylinder_character(
    lambda_partition: Partition, mu_partition: Partition
) -> Character:
    character: Character = {}
    lambda_boxes = partition_boxes(lambda_partition)
    mu_boxes = partition_boxes(mu_partition)
    for j, k in lambda_boxes:
        add_terms(character, (0, j, k), YZ_TERMS, 1)
    for i, k in mu_boxes:
        add_terms(character, (i, 0, k), XZ_TERMS, 1)
    for j, lambda_row in lambda_boxes:
        for i, mu_row in mu_boxes:
            if lambda_row == mu_row:
                add_terms(
                    character,
                    (i, j, lambda_row),
                    CUBE_TERMS,
                    -1,
                )
    return character


def quotient_character(configuration: frozenset[Box]) -> Character:
    character: Character = {}
    for _, i, j, k in configuration:
        add_terms(character, (i, j, k), CUBE_TERMS, 1)
    return character


def merge_characters(left: Character, right: Character) -> Character:
    merged = dict(left)
    for weight, coefficient in right.items():
        add_character(merged, weight, coefficient)
    return merged


def descendent_coefficient(
    character: Character,
    degree: int,
    weights: tuple[Fraction, Fraction, Fraction] = (
        Fraction(1),
        Fraction(2),
        Fraction(-3),
    ),
) -> Fraction:
    total = Fraction(0)
    denominator = factorial(degree)
    for (i, j, k), multiplicity in character.items():
        weight = weights[0] * i + weights[1] * j + weights[2] * k
        total += ((-weight) ** degree) * multiplicity / denominator
    return total


def enumerate_relative_partition_pair_series(
    lambda_partition: Partition,
    mu_partition: Partition,
    maximum_length: int,
    degree: int = 3,
    weights: tuple[Fraction, Fraction, Fraction] = (
        Fraction(1),
        Fraction(2),
        Fraction(-3),
    ),
) -> dict[str, object]:
    """Return exact coefficients of the series P and N through the cutoff."""

    configurations = enumerate_configurations(
        lambda_partition, mu_partition, maximum_length
    )
    no_descendent = [Fraction(0)] * (maximum_length + 1)
    descendent = [Fraction(0)] * (maximum_length + 1)
    base_character = cylinder_character(lambda_partition, mu_partition)
    for configuration in configurations:
        length = len(configuration)
        if length > maximum_length:
            continue
        sign = 1 if length % 2 == 0 else -1
        no_descendent[length] += sign
        character = merge_characters(
            base_character, quotient_character(configuration)
        )
        descendent[length] += sign * descendent_coefficient(
            character, degree, weights
        )
    return {
        "P_coefficients": tuple(no_descendent),
        "N_coefficients": tuple(descendent),
        "configuration_count": len(configurations),
    }


def solve_linear_system(
    matrix: list[list[Fraction]], vector: list[Fraction]
) -> tuple[Fraction, ...] | None:
    dimension = len(matrix)
    augmented = [
        [Fraction(value) for value in row] + [Fraction(vector[index])]
        for index, row in enumerate(matrix)
    ]
    rank = 0
    for column in range(dimension):
        if rank >= dimension:
            break
        pivot = next(
            (
                row
                for row in range(rank, dimension)
                if augmented[row][column] != 0
            ),
            None,
        )
        if pivot is None:
            continue
        augmented[rank], augmented[pivot] = (
            augmented[pivot],
            augmented[rank],
        )
        pivot_value = augmented[rank][column]
        for index in range(column, dimension + 1):
            augmented[rank][index] /= pivot_value
        for row in range(dimension):
            if row == rank or augmented[row][column] == 0:
                continue
            factor = augmented[row][column]
            for index in range(column, dimension + 1):
                augmented[row][index] -= factor * augmented[rank][index]
        rank += 1
    if rank < dimension:
        return None
    return tuple(row[dimension] for row in augmented)


def reconstruct_rational_series(
    coefficients: tuple[Fraction, ...],
    maximum_denominator_degree: int = 18,
    maximum_numerator_degree: int = 12,
    guard: int = 5,
) -> RationalFunction:
    """Reconstruct and guard-check the smallest rational generating series."""

    count = len(coefficients)
    for total_degree in range(
        1, maximum_denominator_degree + maximum_numerator_degree + 2
    ):
        for denominator_degree in range(
            min(maximum_denominator_degree, total_degree - 1) + 1
        ):
            numerator_degree = total_degree - denominator_degree - 1
            if not 0 <= numerator_degree <= maximum_numerator_degree:
                continue
            unknowns = denominator_degree + numerator_degree + 1
            if count < unknowns + guard:
                continue

            matrix: list[list[Fraction]] = []
            vector: list[Fraction] = []
            for coefficient_index in range(unknowns):
                row = [
                    (
                        coefficients[coefficient_index - offset]
                        if coefficient_index - offset >= 0
                        else Fraction(0)
                    )
                    for offset in range(1, denominator_degree + 1)
                ]
                row.extend(
                    Fraction(-1 if coefficient_index == degree else 0)
                    for degree in range(numerator_degree + 1)
                )
                matrix.append(row)
                vector.append(-coefficients[coefficient_index])

            solution = solve_linear_system(matrix, vector)
            if solution is None:
                continue
            valid = True
            for index in range(count):
                left = coefficients[index]
                for offset in range(1, denominator_degree + 1):
                    if index - offset >= 0:
                        left += (
                            solution[offset - 1]
                            * coefficients[index - offset]
                        )
                right = (
                    solution[denominator_degree + index]
                    if index <= numerator_degree
                    else Fraction(0)
                )
                if left != right:
                    valid = False
                    break
            if valid:
                return RationalFunction(
                    solution[denominator_degree:],
                    (Fraction(1),)
                    + solution[:denominator_degree],
                )
    raise ValueError(
        f"could not reconstruct a rational series from {count} coefficients"
    )
