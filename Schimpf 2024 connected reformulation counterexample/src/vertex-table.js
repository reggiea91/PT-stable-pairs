import {
  divideRationals,
  rationalEquals,
  rationalToJSON,
} from "./exact.js";
import {
  enumerateRelativePartitionPairSeries,
  reconstructRationalSeries,
} from "./box-enumerator.js";

export const VERTEX_RELATIVE_PARTITION_PAIRS = Object.freeze([
  Object.freeze({ lambda: Object.freeze([]), mu: Object.freeze([]) }),
  Object.freeze({ lambda: Object.freeze([1]), mu: Object.freeze([]) }),
  Object.freeze({ lambda: Object.freeze([]), mu: Object.freeze([1]) }),
  Object.freeze({ lambda: Object.freeze([1]), mu: Object.freeze([1]) }),
  Object.freeze({ lambda: Object.freeze([2]), mu: Object.freeze([]) }),
  Object.freeze({ lambda: Object.freeze([1, 1]), mu: Object.freeze([]) }),
  Object.freeze({ lambda: Object.freeze([2]), mu: Object.freeze([1]) }),
  Object.freeze({ lambda: Object.freeze([1, 1]), mu: Object.freeze([1]) }),
]);

export function relativePartitionPairKey(lambda, mu) {
  return JSON.stringify([lambda, mu]);
}

export function relativePartitionPairLabel(lambda, mu) {
  const partition = (value) =>
    value.length === 0 ? "empty" : `(${value.join(",")})`;
  return `${partition(lambda)}|${partition(mu)}`;
}

export function deriveVertexTable({
  maximumLength = 18,
  degree = 3,
  weights,
} = {}) {
  const rows = [];
  for (const { lambda, mu } of VERTEX_RELATIVE_PARTITION_PAIRS) {
    const series = enumerateRelativePartitionPairSeries(
      lambda,
      mu,
      maximumLength,
      degree,
      weights
    );
    const P = reconstructRationalSeries(series.noDescendent, {
      maximumDenominatorDegree: 14,
      maximumNumeratorDegree: 10,
      guard: 5,
    });
    const N = reconstructRationalSeries(series.descendent, {
      maximumDenominatorDegree: 18,
      maximumNumeratorDegree: 12,
      guard: 5,
    });
    const h = divideRationals(N, P);
    rows.push({
      key: relativePartitionPairKey(lambda, mu),
      label: relativePartitionPairLabel(lambda, mu),
      lambda: [...lambda],
      mu: [...mu],
      configurationCount: series.configurationCount,
      coefficients: {
        P: series.noDescendent,
        N: series.descendent,
      },
      P,
      N,
      h,
    });
  }
  return {
    maximumLength,
    degree,
    weights: weights || [1, 2, -3],
    rows,
    byKey: new Map(rows.map((row) => [row.key, row])),
  };
}

export function assertVertexTablesEqual(left, right) {
  for (const pair of VERTEX_RELATIVE_PARTITION_PAIRS) {
    const key = relativePartitionPairKey(pair.lambda, pair.mu);
    const leftRow = left.byKey.get(key);
    const rightRow = right.byKey.get(key);
    if (!leftRow || !rightRow) {
      throw new Error(`missing relative-partition pair ${key}`);
    }
    for (const field of ["P", "N", "h"]) {
      if (!rationalEquals(leftRow[field], rightRow[field])) {
        throw new Error(
          `${field} is not stable for ${relativePartitionPairLabel(
            pair.lambda,
            pair.mu
          )}`
        );
      }
    }
  }
}

export function vertexTableToJSON(table, { includeCoefficients = false } = {}) {
  return {
    maximumLength: table.maximumLength,
    degree: table.degree,
    weights: table.weights,
    rows: table.rows.map((row) => ({
      key: row.key,
      label: row.label,
      lambda: row.lambda,
      mu: row.mu,
      configurationCount: row.configurationCount,
      ...(includeCoefficients
        ? {
            enumeratedCoefficients: {
              P: row.coefficients.P.map(String),
              N: row.coefficients.N.map(String),
            },
          }
        : {}),
      P: rationalToJSON(row.P),
      N: rationalToJSON(row.N),
      h: rationalToJSON(row.h),
    })),
  };
}
