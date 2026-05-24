# Quadratic Sieve (QS) - Integer Factorization

Implementation of the Quadratic Sieve algorithm for integer factorization, following the mini-course requirements.

## Overview

This project implements the **Quadratic Sieve (QS)** method for factoring composite integers. The implementation is educational, focusing on clarity and correctness rather than maximum performance.

## Algorithm Steps

The implementation follows the standard QS algorithm:

### 1. Parameter Estimation (`estimate_parameters_task1`)
Computes QS parameters based on L(N) = exp(√(ln N × ln ln N)):
- **theoretical_B**: Factor base bound B = exp(1/√2 × √(ln N × ln ln N))
- **used_B**: Practical B value (larger for better success rate)
- **required_relations**: Number of B-smooth relations needed ≈ |factor_base| + 25
- **sieve_length**: Length of sieve interval ≈ max(20000, B × 35)

### 2. Factor Base Construction (`build_factor_base`)
Builds the factor base: all primes p ≤ B where N is a quadratic residue mod p.
- Uses Tonelli-Shanks algorithm to find square roots mod p
- Returns roots for each prime: solutions to t² ≡ N (mod p)

### 3. Sieving for B-smooth Relations (`find_b_smooth_relations_task2`)
Finds B-smooth values of f(t) = t² - N:
- Uses **logarithmic sieving** with NumPy for speed
- Subtracts log(p) at positions where f(t) ≡ 0 (mod p)
- Candidates with small log values are verified by exact integer division
- Only truly B-smooth values are kept as relations

### 4. Linear Algebra over GF(2) (`find_dependencies_task3`)
Solves for dependencies in the exponent matrix:
- Converts exponent vectors to parity bits (mod 2)
- Uses incremental Gaussian elimination over GF(2)
- Finds subsets of relations whose product is a perfect square

### 5. GCD Computation (`try_dependency_task4`)
For each dependency, constructs X² ≡ Y² (mod N):
- **X** = product of t values from selected relations
- **Y** = square root of product of f(t) values
- Computes gcd(X - Y, N) and gcd(X + Y, N)
- Returns non-trivial factor if found

## Installation

```bash
# Create virtual environment (optional)
python -m venv .venv
source .venv/bin/activate  # Linux/macOS
# .venv\Scripts\activate   # Windows

# Install dependencies
pip install -r requirements.txt
```

**Requirements:**
- Python 3.8+
- NumPy ≥ 1.24

## Usage

### Basic Usage

Factor a composite number:
```bash
python main.py 1022117
```

Output:
```
[input] N=1022117
[input] bit_length=20
[parameters] bits=20 L(N)≈4.155e+02 theoretical_B=71 used_B=426
[factor base] first primes=[2, 17, 31, 47, 59, 61, 71, ...]
[sieve] round=1 interval=[0, 20000) candidates=3506 new_smooth=231
[linear algebra] solving GF(2) matrix with 231 rows and 35 columns
[congruence] X=1011, Y=2
[gcd] gcd(X-Y, N)=1009
[result] factors=[1009, 1013]
[time] 0.046s
```

### Find One Factor Only

```bash
python main.py 11021 --one-factor
```

Output:
```
[result] non_trivial_factor=103
[result] other_factor=107
[check] 103 * 107 = 11021
```

### Specify Factor Base Bound

```bash
python main.py 1022117 --B 200 --one-factor
```

### Save Intermediate Results

```bash
python main.py 1022117 --save-json trace.json
```

The JSON file contains:
- All parameters (B, L(N), sieve length, etc.)
- Complete factor base
- All B-smooth relations found
- Exponent matrix in binary form
- Congruence X² ≡ Y² (mod N)
- GCD results

### Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `N` | Composite integer to factor | Required |
| `--B` | Factor base bound | Auto-estimated |
| `--max-rounds` | Maximum sieve rounds before increasing B | 8 |
| `--one-factor` | Return one factor instead of full factorization | False |
| `--save-json` | Save intermediate results to JSON file | None |

## Running Tests

```bash
python test_qs.py
```

Tests include:
- Small semiprime: 103 × 107 = 11021
- Medium semiprime: 1009 × 1013 = 1022117

## Project Structure

```
quadratic-sieve/
├── main.py                      # CLI entry point
├── qs_course/
│   ├── __init__.py
│   └── quadratic_sieve.py       # Core QS implementation
├── test_qs.py                   # Test cases
├── requirements.txt             # Dependencies
├── .gitignore
└── README.md
```

## Key Functions

### Core Algorithm Functions

- **`estimate_parameters_task1(N)`** - Computes L(N), theoretical B, and other parameters
- **`build_factor_base(N, B)`** - Constructs factor base with quadratic residues
- **`find_b_smooth_relations_task2(...)`** - Sieves for B-smooth values using logarithms
- **`find_dependencies_task3(row_bits)`** - Gaussian elimination over GF(2)
- **`try_dependency_task4(...)`** - Builds X² ≡ Y² (mod N) and computes GCD

### Main API

- **`quadratic_sieve(N, B=None, ...)`** - Returns one non-trivial factor
- **`factor_integer(N, ...)`** - Recursively factors N completely
