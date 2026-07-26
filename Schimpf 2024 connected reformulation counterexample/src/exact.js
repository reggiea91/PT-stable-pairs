/**
 * Exact arithmetic for rational functions in p with coefficients in Q.
 *
 * Polynomials are coefficient arrays in ascending order:
 *   [a0, a1, ...] = a0 + a1*p + ...
 */

export class Fraction {
  constructor(numerator, denominator = 1n) {
    let n = typeof numerator === "bigint" ? numerator : BigInt(numerator);
    let d = typeof denominator === "bigint" ? denominator : BigInt(denominator);
    if (d === 0n) throw new Error("zero denominator");
    if (d < 0n) {
      n = -n;
      d = -d;
    }
    const divisor = Fraction.gcd(n < 0n ? -n : n, d);
    this.n = n / divisor;
    this.d = d / divisor;
    Object.freeze(this);
  }

  static gcd(a, b) {
    while (b !== 0n) [a, b] = [b, a % b];
    return a || 1n;
  }

  static of(value) {
    return value instanceof Fraction ? value : new Fraction(value);
  }

  add(other) {
    const rhs = Fraction.of(other);
    return new Fraction(this.n * rhs.d + rhs.n * this.d, this.d * rhs.d);
  }

  sub(other) {
    return this.add(Fraction.of(other).neg());
  }

  mul(other) {
    const rhs = Fraction.of(other);
    return new Fraction(this.n * rhs.n, this.d * rhs.d);
  }

  div(other) {
    const rhs = Fraction.of(other);
    if (rhs.n === 0n) throw new Error("division by zero");
    return new Fraction(this.n * rhs.d, this.d * rhs.n);
  }

  neg() {
    return new Fraction(-this.n, this.d);
  }

  eq(other) {
    const rhs = Fraction.of(other);
    return this.n === rhs.n && this.d === rhs.d;
  }

  isZero() {
    return this.n === 0n;
  }

  toJSON() {
    return this.toString();
  }

  toString() {
    return this.d === 1n ? `${this.n}` : `${this.n}/${this.d}`;
  }
}

export const q = (numerator, denominator = 1) =>
  new Fraction(BigInt(numerator), BigInt(denominator));

export function polynomial(values) {
  return trimPolynomial(values.map(Fraction.of));
}

export function trimPolynomial(values) {
  let length = values.length;
  while (length > 1 && values[length - 1].isZero()) length -= 1;
  return values.slice(0, Math.max(1, length));
}

export function polynomialIsZero(values) {
  return trimPolynomial(values).every((coefficient) => coefficient.isZero());
}

export function addPolynomials(left, right) {
  const length = Math.max(left.length, right.length);
  return trimPolynomial(
    Array.from({ length }, (_, index) =>
      (left[index] || q(0)).add(right[index] || q(0))
    )
  );
}

export function negatePolynomial(values) {
  return values.map((coefficient) => coefficient.neg());
}

export function subtractPolynomials(left, right) {
  return addPolynomials(left, negatePolynomial(right));
}

export function scalePolynomial(values, scalar) {
  const factor = Fraction.of(scalar);
  return trimPolynomial(values.map((coefficient) => coefficient.mul(factor)));
}

export function multiplyPolynomials(left, right) {
  const result = Array.from(
    { length: left.length + right.length - 1 },
    () => q(0)
  );
  for (let i = 0; i < left.length; i += 1) {
    for (let j = 0; j < right.length; j += 1) {
      result[i + j] = result[i + j].add(left[i].mul(right[j]));
    }
  }
  return trimPolynomial(result);
}

export function dividePolynomials(dividend, divisor) {
  const denominator = trimPolynomial(divisor);
  if (polynomialIsZero(denominator)) throw new Error("polynomial division by zero");
  let remainder = trimPolynomial(dividend);
  const quotient = Array.from(
    { length: Math.max(1, remainder.length - denominator.length + 1) },
    () => q(0)
  );
  while (
    !polynomialIsZero(remainder) &&
    remainder.length >= denominator.length
  ) {
    const degree = remainder.length - denominator.length;
    const factor = remainder[remainder.length - 1].div(
      denominator[denominator.length - 1]
    );
    quotient[degree] = quotient[degree].add(factor);
    const shifted = Array.from({ length: degree }, () => q(0)).concat(
      scalePolynomial(denominator, factor)
    );
    remainder = subtractPolynomials(remainder, shifted);
  }
  return {
    quotient: trimPolynomial(quotient),
    remainder: trimPolynomial(remainder),
  };
}

export function polynomialGcd(left, right) {
  let a = trimPolynomial(left);
  let b = trimPolynomial(right);
  while (!polynomialIsZero(b)) {
    const { remainder } = dividePolynomials(a, b);
    a = b;
    b = remainder;
  }
  if (polynomialIsZero(a)) return polynomial([1]);
  return scalePolynomial(a, q(1).div(a[a.length - 1]));
}

export function evaluatePolynomial(values, point) {
  const x = Fraction.of(point);
  let result = q(0);
  for (let index = values.length - 1; index >= 0; index -= 1) {
    result = result.mul(x).add(values[index]);
  }
  return result;
}

export function rational(numerator, denominator = polynomial([1])) {
  let top = trimPolynomial(numerator);
  let bottom = trimPolynomial(denominator);
  if (polynomialIsZero(bottom)) throw new Error("zero rational-function denominator");
  if (polynomialIsZero(top)) {
    return Object.freeze({ numerator: polynomial([0]), denominator: polynomial([1]) });
  }
  const divisor = polynomialGcd(top, bottom);
  top = dividePolynomials(top, divisor).quotient;
  bottom = dividePolynomials(bottom, divisor).quotient;
  const normalization = bottom[0].isZero()
    ? bottom[bottom.length - 1]
    : bottom[0];
  top = scalePolynomial(top, q(1).div(normalization));
  bottom = scalePolynomial(bottom, q(1).div(normalization));
  return Object.freeze({ numerator: top, denominator: bottom });
}

export function addRationals(left, right) {
  return rational(
    addPolynomials(
      multiplyPolynomials(left.numerator, right.denominator),
      multiplyPolynomials(right.numerator, left.denominator)
    ),
    multiplyPolynomials(left.denominator, right.denominator)
  );
}

export function subtractRationals(left, right) {
  return rational(
    subtractPolynomials(
      multiplyPolynomials(left.numerator, right.denominator),
      multiplyPolynomials(right.numerator, left.denominator)
    ),
    multiplyPolynomials(left.denominator, right.denominator)
  );
}

export function multiplyRationals(left, right) {
  return rational(
    multiplyPolynomials(left.numerator, right.numerator),
    multiplyPolynomials(left.denominator, right.denominator)
  );
}

export function divideRationals(left, right) {
  if (polynomialIsZero(right.numerator)) throw new Error("rational division by zero");
  return rational(
    multiplyPolynomials(left.numerator, right.denominator),
    multiplyPolynomials(left.denominator, right.numerator)
  );
}

export function scaleRational(value, scalar) {
  return rational(scalePolynomial(value.numerator, scalar), value.denominator);
}

export function rationalEquals(left, right) {
  return polynomialIsZero(
    subtractPolynomials(
      multiplyPolynomials(left.numerator, right.denominator),
      multiplyPolynomials(right.numerator, left.denominator)
    )
  );
}

export function shiftedAtOne(values) {
  const result = Array.from({ length: values.length }, () => q(0));
  for (let degree = 0; degree < values.length; degree += 1) {
    let binomial = 1n;
    for (let index = 0; index <= degree; index += 1) {
      if (index > 0) {
        binomial =
          (binomial * BigInt(degree - index + 1)) / BigInt(index);
      }
      result[index] = result[index].add(
        values[degree].mul(new Fraction(binomial))
      );
    }
  }
  return trimPolynomial(result);
}

export function leadingPolynomialAtOne(values) {
  const shifted = shiftedAtOne(values);
  for (let order = 0; order < shifted.length; order += 1) {
    if (!shifted[order].isZero()) return { order, coefficient: shifted[order] };
  }
  return { order: Number.POSITIVE_INFINITY, coefficient: q(0) };
}

export function leadingRationalAtOne(value) {
  const top = leadingPolynomialAtOne(value.numerator);
  const bottom = leadingPolynomialAtOne(value.denominator);
  return {
    order: top.order - bottom.order,
    coefficient: top.coefficient.div(bottom.coefficient),
  };
}

export function laurentAtOne(value, minimumPower, maximumPower) {
  const shiftedNumerator = shiftedAtOne(value.numerator);
  const shiftedDenominator = shiftedAtOne(value.denominator);
  const numeratorLead = leadingPolynomialAtOne(value.numerator);
  const denominatorLead = leadingPolynomialAtOne(value.denominator);
  if (numeratorLead.order === Number.POSITIVE_INFINITY) {
    return new Map(
      Array.from(
        { length: maximumPower - minimumPower + 1 },
        (_, index) => [minimumPower + index, q(0)]
      )
    );
  }
  const start = numeratorLead.order - denominatorLead.order;
  const numerator = shiftedNumerator.slice(numeratorLead.order);
  const denominator = shiftedDenominator.slice(denominatorLead.order);
  const needed = Math.max(0, maximumPower - start);
  const coefficients = [];
  for (let index = 0; index <= needed; index += 1) {
    let rhs = numerator[index] || q(0);
    for (let offset = 1; offset <= index; offset += 1) {
      rhs = rhs.sub(
        (denominator[offset] || q(0)).mul(coefficients[index - offset])
      );
    }
    coefficients[index] = rhs.div(denominator[0]);
  }
  const output = new Map();
  for (let power = minimumPower; power <= maximumPower; power += 1) {
    const index = power - start;
    output.set(
      power,
      index >= 0 && index < coefficients.length ? coefficients[index] : q(0)
    );
  }
  return output;
}

export function polynomialToString(values, variable = "p") {
  const terms = [];
  for (let degree = 0; degree < values.length; degree += 1) {
    const coefficient = values[degree];
    if (coefficient.isZero()) continue;
    const negative = coefficient.n < 0n;
    const magnitude = new Fraction(
      negative ? -coefficient.n : coefficient.n,
      coefficient.d
    ).toString();
    const variablePart =
      degree === 1 ? variable : `${variable}^${degree}`;
    const body =
      degree === 0
        ? magnitude
        : magnitude === "1"
          ? variablePart
          : `${magnitude}*${variablePart}`;
    terms.push({ sign: negative ? "-" : "+", body });
  }
  if (terms.length === 0) return "0";
  const [first, ...rest] = terms;
  let output = first.sign === "-" ? `-${first.body}` : first.body;
  for (const term of rest) output += ` ${term.sign} ${term.body}`;
  return output;
}

export function rationalToString(value, variable = "p") {
  const numerator = polynomialToString(value.numerator, variable);
  const denominator = polynomialToString(value.denominator, variable);
  return denominator === "1" ? numerator : `(${numerator})/(${denominator})`;
}

export function rationalToJSON(value) {
  return {
    numerator: value.numerator.map(String),
    denominator: value.denominator.map(String),
    expression: rationalToString(value),
  };
}
