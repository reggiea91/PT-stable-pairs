import {
  addRationals,
  laurentAtOne,
  leadingRationalAtOne,
  multiplyPolynomials,
  multiplyRationals,
  polynomial,
  q,
  rational,
  rationalEquals,
  rationalToJSON,
  scaleRational,
  subtractRationals,
} from "./exact.js";
import { relativePartitionPairKey } from "./vertex-table.js";

function polynomialPower(base, exponent) {
  let result = polynomial([1]);
  for (let index = 0; index < exponent; index += 1) {
    result = multiplyPolynomials(result, base);
  }
  return result;
}

/**
 * Toric-chain no-descendent gluing weights A_{lambda|mu}.
 *
 * These are the geometric edge/gluing factors. The vertex ratios h that they
 * multiply are derived independently by the box enumerator.
 */
export function toricChainWeightTable() {
  const pPlusOne = polynomial([1, 1]);
  const pMinusOne = polynomial([-1, 1]);
  const pPlusOne2 = polynomialPower(pPlusOne, 2);
  const pPlusOne4 = polynomialPower(pPlusOne, 4);
  const pPlusOne6 = polynomialPower(pPlusOne, 6);
  const pMinusOne2 = polynomialPower(pMinusOne, 2);
  return new Map([
    [relativePartitionPairKey([], []), rational(polynomial([1]))],
    [relativePartitionPairKey([1], []), rational(polynomial([0, 1]), pPlusOne2)],
    [relativePartitionPairKey([], [1]), rational(polynomial([0, 1]), pPlusOne2)],
    [
      relativePartitionPairKey([1], [1]),
      rational(polynomial([0, -1, -1, -1]), pPlusOne4),
    ],
    [
      relativePartitionPairKey([2], []),
      rational(
        polynomial([0, 0, 0, -1]),
        multiplyPolynomials(pMinusOne2, pPlusOne4)
      ),
    ],
    [
      relativePartitionPairKey([1, 1], []),
      rational(
        polynomial([0, 0, 0, -1]),
        multiplyPolynomials(pMinusOne2, pPlusOne4)
      ),
    ],
    [
      relativePartitionPairKey([2], [1]),
      rational(
        polynomial([0, 0, -1, -1, 0, 1]),
        multiplyPolynomials(pMinusOne2, pPlusOne6)
      ),
    ],
    [
      relativePartitionPairKey([1, 1], [1]),
      rational(
        polynomial([0, 0, 0, 1, 0, -1, -1]),
        multiplyPolynomials(pMinusOne2, pPlusOne6)
      ),
    ],
  ]);
}

/**
 * Feed the enumerated vertex table directly into the connected quotient.
 */
export function computeConnectedQuotient(vertexTable) {
  const A = toricChainWeightTable();
  const B = new Map(
    vertexTable.rows.map((row) => {
      const weight = A.get(row.key);
      if (!weight) throw new Error(`missing toric-chain weight for ${row.key}`);
      return [row.key, multiplyRationals(weight, row.h)];
    })
  );
  const getA = (lambda, mu) => A.get(relativePartitionPairKey(lambda, mu));
  const getB = (lambda, mu) => B.get(relativePartitionPairKey(lambda, mu));

  const Z10 = getA([1], []);
  const Z01 = getA([], [1]);
  const Z11 = getA([1], [1]);
  const Z20 = addRationals(getA([2], []), getA([1, 1], []));

  const F10 = getB([1], []);
  const F01 = getB([], [1]);
  const F11 = getB([1], [1]);
  const F20 = addRationals(getB([2], []), getB([1, 1], []));
  const F21 = addRationals(getB([2], [1]), getB([1, 1], [1]));

  const C11 = subtractRationals(
    subtractRationals(F11, multiplyRationals(Z10, F01)),
    multiplyRationals(Z01, F10)
  );

  const terms = [
    { label: "F_{2,1}", value: F21 },
    {
      label: "-Z_{0,1}F_{2,0}",
      value: scaleRational(multiplyRationals(Z01, F20), q(-1)),
    },
    {
      label: "-Z_{1,0}F_{1,1}",
      value: scaleRational(multiplyRationals(Z10, F11), q(-1)),
    },
    {
      label: "-Z_{2,0}F_{0,1}",
      value: scaleRational(multiplyRationals(Z20, F01), q(-1)),
    },
    {
      label: "-Z_{1,1}F_{1,0}",
      value: scaleRational(multiplyRationals(Z11, F10), q(-1)),
    },
    {
      label: "2Z_{1,0}Z_{0,1}F_{1,0}",
      value: scaleRational(
        multiplyRationals(multiplyRationals(Z10, Z01), F10),
        q(2)
      ),
    },
    {
      label: "Z_{1,0}^2F_{0,1}",
      value: multiplyRationals(multiplyRationals(Z10, Z10), F01),
    },
  ];

  const C21 = terms.reduce(
    (total, term) => addRationals(total, term.value),
    rational(polynomial([0]))
  );
  const expectedC21 = rational(
    polynomial([0, 0, 18]),
    polynomial([-1, -2, 0, 2, 1])
  );
  if (!rationalEquals(C21, expectedC21)) {
    throw new Error(
      "connected quotient does not simplify to 18*p^2/((p-1)*(p+1)^3)"
    );
  }

  const laurentTerms = terms.map((term) => ({
    ...term,
    laurent: laurentAtOne(term.value, -3, 0),
  }));
  const laurent = laurentAtOne(C21, -3, 0);

  return {
    A,
    B,
    aggregates: { Z10, Z01, Z11, Z20, F10, F01, F11, F20, F21 },
    C11,
    C21,
    expectedC21,
    C11Leading: leadingRationalAtOne(C11),
    C21Leading: leadingRationalAtOne(C21),
    terms: laurentTerms,
    laurent,
  };
}

export function connectedQuotientToJSON(result) {
  const serializeMap = (map) =>
    Object.fromEntries([...map].map(([key, value]) => [key, rationalToJSON(value)]));
  const serializeLaurent = (map) =>
    Object.fromEntries([...map].map(([power, coefficient]) => [power, String(coefficient)]));
  return {
    toricChainWeights: serializeMap(result.A),
    vertexContributions: serializeMap(result.B),
    aggregates: Object.fromEntries(
      Object.entries(result.aggregates).map(([key, value]) => [
        key,
        rationalToJSON(value),
      ])
    ),
    C11: rationalToJSON(result.C11),
    C21: rationalToJSON(result.C21),
    C21Expected: rationalToJSON(result.expectedC21),
    C11LeadingAtOne: {
      order: Number.isFinite(result.C11Leading.order)
        ? result.C11Leading.order
        : "Infinity",
      coefficient: String(result.C11Leading.coefficient),
    },
    C21LeadingAtOne: {
      order: result.C21Leading.order,
      coefficient: String(result.C21Leading.coefficient),
    },
    laurentAtOne: serializeLaurent(result.laurent),
    laurentTerms: result.terms.map((term) => ({
      label: term.label,
      value: rationalToJSON(term.value),
      laurentAtOne: serializeLaurent(term.laurent),
    })),
  };
}
