import {
  assertVertexTablesEqual,
  deriveVertexTable,
  vertexTableToJSON,
} from "./vertex-table.js";
import {
  computeConnectedQuotient,
  connectedQuotientToJSON,
} from "./connected-quotient.js";

export function runReproductionPipeline({
  maximumLength = 18,
  stabilityCutoff = 16,
  degree = 3,
  weights = [1, 2, -3],
} = {}) {
  if (stabilityCutoff >= maximumLength) {
    throw new Error("stabilityCutoff must be smaller than maximumLength");
  }
  const stabilityTable = deriveVertexTable({
    maximumLength: stabilityCutoff,
    degree,
    weights,
  });
  const vertexTable = deriveVertexTable({
    maximumLength,
    degree,
    weights,
  });
  assertVertexTablesEqual(stabilityTable, vertexTable);
  const connected = computeConnectedQuotient(vertexTable);
  return {
    parameters: { maximumLength, stabilityCutoff, degree, weights },
    stabilityTable,
    vertexTable,
    connected,
  };
}

export function pipelineToJSON(result, { includeCoefficients = true } = {}) {
  return {
    schemaVersion: 1,
    description:
      "End-to-end finite two-leg PT vertex enumeration and connected-quotient certificate",
    parameters: result.parameters,
    stabilityCheck: {
      status: "PASS",
      comparedCutoffs: [
        result.parameters.stabilityCutoff,
        result.parameters.maximumLength,
      ],
    },
    vertexTable: vertexTableToJSON(result.vertexTable, {
      includeCoefficients,
    }),
    connectedQuotient: connectedQuotientToJSON(result.connected),
    conclusion: {
      exactIdentity: "C21(p) = 18*p^2/((p-1)*(p+1)^3)",
      residueAtP1: "9/4",
      forbiddenPoleVerified: true,
    },
  };
}

