# Mini-course: Factorization by Quadratic Sieve

This project implements integer factorization by the **Quadratic Sieve (QS)** method.

## What is implemented

The code follows the required mini-course steps:

1. Estimate QS parameters:
   - factor-base bound `B`,
   - required number of B-smooth relations,
   - expected sieve interval length.
2. Build a factor base of primes `p <= B` for which `N` is a quadratic residue modulo `p`.
3. Search for B-smooth values of

   ```text
   Q(x) = x^2 - N
   ```

   using logarithmic sieving.
4. Build exponent vectors modulo 2.
5. Solve linear dependencies over `GF(2)`.
6. Construct

   ```text
   X^2 ≡ Y^2 (mod N)
   ```

7. Compute

   ```text
   gcd(X - Y, N)
   gcd(X + Y, N)
   ```

   to obtain a non-trivial divisor of `N`.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate      # Linux/macOS
# .venv\Scripts\activate      # Windows PowerShell
pip install -r requirements.txt
```

## Run

```bash
python main.py 635710871264254320370996442119 --one-factor
```

Or full recursive factorization:

```bash
python main.py 635710871264254320370996442119
```

Save intermediate results for defense:

```bash
python main.py 635710871264254320370996442119 --one-factor --save-json trace.json
```

## Useful options

```bash
python main.py N --B 30000 --max-rounds 10 --one-factor --save-json trace.json
```

- `--B`: manually set factor-base bound.
- `--max-rounds`: number of sieve intervals before increasing `B`.
- `--one-factor`: return one non-trivial factor using QS.
- `--save-json`: save intermediate data: parameters, factor base, smooth relations, matrix rows, congruence, GCD results.

## Notes for oral defense

Important functions:

- `estimate_parameters()` — computes `L(N)`, theoretical `B`, used `B`, relation target, sieve length.
- `build_factor_base()` — constructs the factor base.
- `sieve_interval()` — finds B-smooth values with logarithmic sieving.
- `factor_over_base()` — exact verification of B-smoothness.
- `find_dependencies()` — Gaussian elimination over `GF(2)` using bitsets.
- `try_dependency()` — creates `X^2 ≡ Y^2 (mod N)` and applies `gcd`.
- `quadratic_sieve()` — main QS algorithm.

## Why NumPy is used

NumPy is used only to speed up the logarithmic sieve. All final smoothness checks and GCD calculations use exact integer arithmetic.
