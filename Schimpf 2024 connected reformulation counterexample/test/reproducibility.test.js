import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

import {
  divideRationals,
  polynomial,
  rational,
  rationalEquals,
} from "../src/exact.js";
import {
  assertVertexTablesEqual,
  deriveVertexTable,
  relativePartitionPairKey,
} from "../src/vertex-table.js";
import { computeConnectedQuotient } from "../src/connected-quotient.js";
import { pipelineToJSON, runReproductionPipeline } from "../src/pipeline.js";

const R = (numerator, denominator = [1]) =>
  rational(polynomial(numerator), polynomial(denominator));

function assertRationalEquals(actual, expected, message) {
  assert.ok(rationalEquals(actual, expected), message);
}

const EXPECTED_VERTEX_ROWS = new Map([
  [
    relativePartitionPairKey([], []),
    {
      P: R([1]),
      N: R([0]),
    },
  ],
  [
    relativePartitionPairKey([1], []),
    {
      P: R([1], [1, 1]),
      N: R([-3, 3], [1, 2, 1]),
    },
  ],
  [
    relativePartitionPairKey([], [1]),
    {
      P: R([1], [1, 1]),
      N: R([-3, 3], [1, 2, 1]),
    },
  ],
  [
    relativePartitionPairKey([1], [1]),
    {
      P: R([1, 1, 1], [1, 2, 1]),
      N: R([0, 6, -6], [1, 3, 3, 1]),
    },
  ],
  [
    relativePartitionPairKey([2], []),
    {
      P: R([1], [1, 1, -1, -1]),
      N: R([6, 6, -24], [1, 1, -2, -2, 1, 1]),
    },
  ],
  [
    relativePartitionPairKey([1, 1], []),
    {
      P: R([1], [1, 1, -1, -1]),
      N: R([-24, 6, 6], [1, 1, -2, -2, 1, 1]),
    },
  ],
  [
    relativePartitionPairKey([2], [1]),
    {
      P: R([1, 1, 0, -1], [1, 2, 0, -2, -1]),
      N: R(
        [15, 21, -27, -30, -12, 21],
        [1, 2, -1, -4, -1, 2, 1]
      ),
    },
  ],
  [
    relativePartitionPairKey([1, 1], [1]),
    {
      P: R([1, 0, -1, -1], [1, 2, 0, -2, -1]),
      N: R(
        [-21, 12, 30, 27, -21, -15],
        [1, 2, -1, -4, -1, 2, 1]
      ),
    },
  ],
]);

test("exact rational functions are reduced canonically", () => {
  const value = R([0, 0, 18], [-1, -2, 0, 2, 1]);
  assert.deepEqual(value.numerator.map(String), ["0", "0", "-18"]);
  assert.deepEqual(value.denominator.map(String), ["1", "2", "0", "-2", "-1"]);
});

test("box enumeration reconstructs every P and N row printed in the paper", () => {
  const table = deriveVertexTable({ maximumLength: 18 });
  for (const [key, expected] of EXPECTED_VERTEX_ROWS) {
    const row = table.byKey.get(key);
    assert.ok(row, `missing relative-partition pair ${key}`);
    assertRationalEquals(row.P, expected.P, `P mismatch for ${key}`);
    assertRationalEquals(row.N, expected.N, `N mismatch for ${key}`);
    assertRationalEquals(
      row.h,
      divideRationals(expected.N, expected.P),
      `h=N/P mismatch for ${key}`
    );
  }
});

test("rational reconstruction is stable between cutoffs 16 and 18", () => {
  const cutoff16 = deriveVertexTable({ maximumLength: 16 });
  const cutoff18 = deriveVertexTable({ maximumLength: 18 });
  assert.doesNotThrow(() => assertVertexTablesEqual(cutoff16, cutoff18));
});

test("enumerated vertex table feeds the connected quotient", () => {
  const table = deriveVertexTable({ maximumLength: 18 });
  const result = computeConnectedQuotient(table);
  assertRationalEquals(result.C11, R([0]), "C11 must vanish");
  assertRationalEquals(
    result.C21,
    R([0, 0, 18], [-1, -2, 0, 2, 1]),
    "C21 exact identity mismatch"
  );
  assert.equal(result.C21Leading.order, -1);
  assert.equal(String(result.C21Leading.coefficient), "9/4");
  assert.equal(String(result.laurent.get(-3)), "0");
  assert.equal(String(result.laurent.get(-2)), "0");
  assert.equal(String(result.laurent.get(-1)), "9/4");
  assert.equal(String(result.laurent.get(0)), "9/8");
});

test("committed JSON certificate matches a fresh end-to-end run", async () => {
  const pipeline = runReproductionPipeline();
  const fresh = pipelineToJSON(pipeline, { includeCoefficients: true });
  const committed = JSON.parse(
    await readFile(
      new URL("../results/certificate.json", import.meta.url),
      "utf8"
    )
  );
  assert.deepEqual(committed, fresh);
});
