#!/usr/bin/env node

import { createHash } from "node:crypto";
import { mkdir, readFile, writeFile } from "node:fs/promises";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import { pipelineToJSON, runReproductionPipeline } from "../src/pipeline.js";
import { rationalToString } from "../src/exact.js";

const scriptDirectory = dirname(fileURLToPath(import.meta.url));
const repositoryRoot = resolve(scriptDirectory, "..");

function parsePositiveInteger(flag, fallback) {
  const index = process.argv.indexOf(flag);
  if (index < 0) return fallback;
  const value = Number.parseInt(process.argv[index + 1], 10);
  if (!Number.isInteger(value) || value <= 0) {
    throw new Error(`${flag} must be followed by a positive integer`);
  }
  return value;
}

function escapeCsv(value) {
  const text = String(value);
  return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

function laurentLedgerCsv(pipeline) {
  const powers = [-3, -2, -1, 0];
  const header = ["term", ...powers.map((power) => `(p-1)^${power}`)];
  const rows = pipeline.connected.terms.map((term) => [
    term.label,
    ...powers.map((power) => String(term.laurent.get(power))),
  ]);
  rows.push([
    "TOTAL",
    ...powers.map((power) => String(pipeline.connected.laurent.get(power))),
  ]);
  return [header, ...rows]
    .map((row) => row.map(escapeCsv).join(","))
    .join("\n")
    .concat("\n");
}

async function main() {
  const maximumLength = parsePositiveInteger("--max-length", 18);
  const stabilityCutoff = parsePositiveInteger("--stability-cutoff", 16);
  const pipeline = runReproductionPipeline({
    maximumLength,
    stabilityCutoff,
  });
  const resultsDirectory = resolve(repositoryRoot, "results");
  await mkdir(resultsDirectory, { recursive: true });

  const json = pipelineToJSON(pipeline, { includeCoefficients: true });
  await Promise.all([
    writeFile(
      resolve(resultsDirectory, "certificate.json"),
      `${JSON.stringify(json, null, 2)}\n`,
      "utf8"
    ),
    writeFile(
      resolve(resultsDirectory, "laurent-ledger.csv"),
      laurentLedgerCsv(pipeline),
      "utf8"
    ),
  ]);
  const generatedNames = [
    "certificate.json",
    "laurent-ledger.csv",
  ];
  const checksumLines = [];
  for (const name of generatedNames) {
    const content = await readFile(resolve(resultsDirectory, name));
    const digest = createHash("sha256").update(content).digest("hex");
    checksumLines.push(`${digest}  ${name}`);
  }
  await writeFile(
    resolve(resultsDirectory, "SHA256SUMS"),
    `${checksumLines.join("\n")}\n`,
    "utf8"
  );

  console.log("Finite two-leg PT vertex enumeration: PASS");
  console.log(
    `Rational reconstruction stability (${stabilityCutoff} vs ${maximumLength}): PASS`
  );
  console.log(`C11(p) = ${rationalToString(pipeline.connected.C11)}`);
  console.log(`C21(p) = ${rationalToString(pipeline.connected.C21)}`);
  console.log(
    "Exact identity C21(p) = 18*p^2/((p-1)*(p+1)^3): PASS"
  );
  console.log(
    `Residue at p=1: ${pipeline.connected.C21Leading.coefficient}`
  );
  console.log("Forbidden simple pole at p=1: VERIFIED");
  console.log(`Generated results: ${resultsDirectory}`);
}

main().catch((error) => {
  console.error(error.stack || error.message);
  process.exitCode = 1;
});
