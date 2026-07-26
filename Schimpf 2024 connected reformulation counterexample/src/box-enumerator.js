import {
  Fraction,
  q,
  rational,
} from "./exact.js";

function partitionBoxes(partition) {
  const boxes = [];
  for (let row = 0; row < partition.length; row += 1) {
    for (let column = 0; column < partition[row]; column += 1) {
      boxes.push([column, row]);
    }
  }
  return boxes;
}

const boxKey = (box) => box.join(",");
const parseBox = (key) => key.split(",").map(Number);

function containsPair(pairs, first, second) {
  return pairs.some(([left, right]) => left === first && right === second);
}

function memberships([i, j, k], lambda, mu) {
  const membershipsAtWeight = [];
  if (containsPair(partitionBoxes(lambda), j, k)) membershipsAtWeight.push(1);
  if (containsPair(partitionBoxes(mu), i, k)) membershipsAtWeight.push(2);
  return membershipsAtWeight;
}

function boxesAtWeight(weight, lambda, mu) {
  const atWeight = memberships(weight, lambda, mu);
  if (weight.some((coordinate) => coordinate < 0)) {
    return atWeight.length === 1 ? [[atWeight[0], ...weight]] : [];
  }
  return atWeight.length === 2 ? [[12, ...weight]] : [];
}

function quotientBoxes(lambda, mu, maximumLength) {
  const boxes = new Map();
  for (const [j, k] of partitionBoxes(lambda)) {
    for (let i = -maximumLength; i < 0; i += 1) {
      boxes.set(boxKey([1, i, j, k]), [1, i, j, k]);
    }
  }
  for (const [i, k] of partitionBoxes(mu)) {
    for (let j = -maximumLength; j < 0; j += 1) {
      boxes.set(boxKey([2, i, j, k]), [2, i, j, k]);
    }
  }
  for (const [j, lambdaRow] of partitionBoxes(lambda)) {
    for (const [i, muRow] of partitionBoxes(mu)) {
      if (lambdaRow === muRow) {
        boxes.set(boxKey([12, i, j, lambdaRow]), [12, i, j, lambdaRow]);
      }
    }
  }
  return [...boxes.values()];
}

function successors(box, lambda, mu) {
  const [, i, j, k] = box;
  const output = [];
  for (const [di, dj, dk] of [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
  ]) {
    output.push(...boxesAtWeight([i + di, j + dj, k + dk], lambda, mu));
  }
  return output;
}

/**
 * Enumerate upward-closed finite quotient-box configurations.
 */
export function enumerateConfigurations(lambda, mu, maximumLength) {
  const boxes = quotientBoxes(lambda, mu, maximumLength);
  const available = new Set(boxes.map(boxKey));
  const successorMap = new Map(
    boxes.map((box) => [
      boxKey(box),
      successors(box, lambda, mu)
        .filter((successor) => available.has(boxKey(successor)))
        .map(boxKey),
    ])
  );
  const ordered = [...boxes].sort(
    (left, right) =>
      right[1] +
        right[2] +
        right[3] -
        (left[1] + left[2] + left[3]) ||
      boxKey(right).localeCompare(boxKey(left))
  );
  const chosen = new Set();
  const configurations = [];

  function visit(index) {
    if (chosen.size > maximumLength) return;
    if (index === ordered.length) {
      configurations.push([...chosen]);
      return;
    }
    const key = boxKey(ordered[index]);
    visit(index + 1);
    if (successorMap.get(key).every((successor) => chosen.has(successor))) {
      chosen.add(key);
      visit(index + 1);
      chosen.delete(key);
    }
  }

  visit(0);
  return configurations;
}

function addCharacter(character, weight, coefficient) {
  const key = weight.join(",");
  const next = (character.get(key) || 0) + coefficient;
  if (next === 0) character.delete(key);
  else character.set(key, next);
}

const cubeTerms = [
  [0, 0, 0, 1],
  [1, 0, 0, -1],
  [0, 1, 0, -1],
  [0, 0, 1, -1],
  [1, 1, 0, 1],
  [1, 0, 1, 1],
  [0, 1, 1, 1],
  [1, 1, 1, -1],
];

const yzTerms = [
  [0, 0, 0, 1],
  [0, 1, 0, -1],
  [0, 0, 1, -1],
  [0, 1, 1, 1],
];

const xzTerms = [
  [0, 0, 0, 1],
  [1, 0, 0, -1],
  [0, 0, 1, -1],
  [1, 0, 1, 1],
];

function addTerms(character, base, terms, coefficient) {
  for (const [di, dj, dk, sign] of terms) {
    addCharacter(
      character,
      [base[0] + di, base[1] + dj, base[2] + dk],
      coefficient * sign
    );
  }
}

function cylinderCharacter(lambda, mu) {
  const character = new Map();
  for (const [j, k] of partitionBoxes(lambda)) {
    addTerms(character, [0, j, k], yzTerms, 1);
  }
  for (const [i, k] of partitionBoxes(mu)) {
    addTerms(character, [i, 0, k], xzTerms, 1);
  }
  for (const [j, lambdaRow] of partitionBoxes(lambda)) {
    for (const [i, muRow] of partitionBoxes(mu)) {
      if (lambdaRow === muRow) {
        addTerms(character, [i, j, lambdaRow], cubeTerms, -1);
      }
    }
  }
  return character;
}

function quotientCharacter(configuration) {
  const character = new Map();
  for (const key of configuration) {
    const [, i, j, k] = parseBox(key);
    addTerms(character, [i, j, k], cubeTerms, 1);
  }
  return character;
}

function mergeCharacters(left, right) {
  const merged = new Map(left);
  for (const [key, coefficient] of right) {
    const next = (merged.get(key) || 0) + coefficient;
    if (next === 0) merged.delete(key);
    else merged.set(key, next);
  }
  return merged;
}

function factorial(number) {
  let result = 1n;
  for (let value = 2n; value <= BigInt(number); value += 1n) result *= value;
  return result;
}

function fractionPower(value, exponent) {
  let result = q(1);
  for (let index = 0; index < exponent; index += 1) result = result.mul(value);
  return result;
}

function descendentCoefficient(
  character,
  degree,
  weights = [q(1), q(2), q(-3)]
) {
  const exactWeights = weights.map(Fraction.of);
  let total = q(0);
  const denominator = new Fraction(factorial(degree));
  for (const [key, multiplicity] of character) {
    const [i, j, k] = key.split(",").map(Number);
    const weight = exactWeights[0]
      .mul(i)
      .add(exactWeights[1].mul(j))
      .add(exactWeights[2].mul(k));
    total = total.add(
      fractionPower(weight.neg(), degree).div(denominator).mul(multiplicity)
    );
  }
  return total;
}

/**
 * Return exact coefficient lists for the no-descendent series P and
 * descendent series N through the requested length cutoff.
 */
export function enumerateRelativePartitionPairSeries(
  lambda,
  mu,
  maximumLength,
  degree = 3,
  weights = [q(1), q(2), q(-3)]
) {
  const configurations = enumerateConfigurations(lambda, mu, maximumLength);
  const noDescendent = Array.from({ length: maximumLength + 1 }, () => q(0));
  const descendent = Array.from({ length: maximumLength + 1 }, () => q(0));
  const baseCharacter = cylinderCharacter(lambda, mu);
  for (const configuration of configurations) {
    const length = configuration.length;
    if (length > maximumLength) continue;
    const sign = length % 2 === 0 ? 1 : -1;
    noDescendent[length] = noDescendent[length].add(sign);
    const character = mergeCharacters(
      baseCharacter,
      quotientCharacter(configuration)
    );
    descendent[length] = descendent[length].add(
      descendentCoefficient(character, degree, weights).mul(sign)
    );
  }
  return {
    noDescendent,
    descendent,
    configurationCount: configurations.length,
  };
}

function solveLinearSystem(matrix, vector) {
  const dimension = matrix.length;
  const augmented = matrix.map((row, index) =>
    row.map(Fraction.of).concat([Fraction.of(vector[index])])
  );
  let rank = 0;
  for (let column = 0; column < dimension && rank < dimension; column += 1) {
    let pivot = -1;
    for (let row = rank; row < dimension; row += 1) {
      if (!augmented[row][column].isZero()) {
        pivot = row;
        break;
      }
    }
    if (pivot < 0) continue;
    [augmented[rank], augmented[pivot]] = [
      augmented[pivot],
      augmented[rank],
    ];
    const pivotValue = augmented[rank][column];
    for (let index = column; index <= dimension; index += 1) {
      augmented[rank][index] = augmented[rank][index].div(pivotValue);
    }
    for (let row = 0; row < dimension; row += 1) {
      if (row === rank || augmented[row][column].isZero()) continue;
      const factor = augmented[row][column];
      for (let index = column; index <= dimension; index += 1) {
        augmented[row][index] = augmented[row][index].sub(
          factor.mul(augmented[rank][index])
        );
      }
    }
    rank += 1;
  }
  if (rank < dimension) return null;
  return augmented.map((row) => row[dimension]);
}

/**
 * Reconstruct a rational generating function from exact initial coefficients.
 * The final `guard` coefficients are verification coefficients not used to fit
 * the smallest candidate.
 */
export function reconstructRationalSeries(
  coefficients,
  { maximumDenominatorDegree = 18, maximumNumeratorDegree = 12, guard = 5 } = {}
) {
  const count = coefficients.length;
  for (
    let totalDegree = 1;
    totalDegree <= maximumDenominatorDegree + maximumNumeratorDegree + 1;
    totalDegree += 1
  ) {
    for (
      let denominatorDegree = 0;
      denominatorDegree <=
      Math.min(maximumDenominatorDegree, totalDegree - 1);
      denominatorDegree += 1
    ) {
      const numeratorDegree = totalDegree - denominatorDegree - 1;
      if (
        numeratorDegree < 0 ||
        numeratorDegree > maximumNumeratorDegree
      ) {
        continue;
      }
      const unknowns = denominatorDegree + numeratorDegree + 1;
      if (count < unknowns + guard) continue;
      const matrix = [];
      const vector = [];
      for (let coefficientIndex = 0; coefficientIndex < unknowns; coefficientIndex += 1) {
        const row = [];
        for (let offset = 1; offset <= denominatorDegree; offset += 1) {
          row.push(
            coefficientIndex - offset >= 0
              ? coefficients[coefficientIndex - offset]
              : q(0)
          );
        }
        for (let degree = 0; degree <= numeratorDegree; degree += 1) {
          row.push(coefficientIndex === degree ? q(-1) : q(0));
        }
        matrix.push(row);
        vector.push(coefficients[coefficientIndex].neg());
      }
      const solution = solveLinearSystem(matrix, vector);
      if (!solution) continue;
      let valid = true;
      for (let index = 0; index < count; index += 1) {
        let left = coefficients[index];
        for (let offset = 1; offset <= denominatorDegree; offset += 1) {
          if (index - offset >= 0) {
            left = left.add(
              solution[offset - 1].mul(coefficients[index - offset])
            );
          }
        }
        const right =
          index <= numeratorDegree
            ? solution[denominatorDegree + index]
            : q(0);
        if (!left.eq(right)) {
          valid = false;
          break;
        }
      }
      if (valid) {
        return rational(
          solution.slice(denominatorDegree),
          [q(1), ...solution.slice(0, denominatorDegree)]
        );
      }
    }
  }
  throw new Error(
    `could not reconstruct a rational series from ${coefficients.length} coefficients`
  );
}
