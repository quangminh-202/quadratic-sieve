from __future__ import annotations

import json
import math
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

try:
    import numpy as np
except Exception as exc:  # pragma: no cover
    np = None
    _NUMPY_IMPORT_ERROR = exc
else:
    _NUMPY_IMPORT_ERROR = None


@dataclass
class QSParameters:
    """Parameters estimated for one input N."""

    n_bits: int
    ln_n: float
    ln_ln_n: float
    L_n: float
    theoretical_B: int
    used_B: int
    factor_base_size: int = 0
    required_relations: int = 0
    sieve_length: int = 0


@dataclass
class Relation:
    """One B-smooth relation f(t) = t^2 - N."""

    t: int
    f_value: int
    exponents: List[int]
    row_bits: int

    def to_light_dict(self) -> dict:
        return {
            "t": self.t,
            "f_value": self.f_value,
            "exponents": self.exponents,
            "row_bits_binary": bin(self.row_bits),
        }


@dataclass
class QSTrace:
    """Data saved for explanation/defense."""

    N: int
    parameters: Optional[QSParameters] = None
    factor_base: List[int] = field(default_factory=list)
    sieve_rounds: List[dict] = field(default_factory=list)
    relations: List[dict] = field(default_factory=list)
    matrix_rows_binary: List[str] = field(default_factory=list)
    dependencies_tried: int = 0
    congruence: Optional[dict] = None
    result: Optional[dict] = None

    def save_json(self, path: str | Path) -> None:
        data = asdict(self)
        if self.parameters is not None:
            data["parameters"] = asdict(self.parameters)
        Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def ceil_sqrt(n: int) -> int:
    r = math.isqrt(n)
    return r if r * r == n else r + 1


def is_square(n: int) -> Tuple[bool, int]:
    if n < 0:
        return False, 0
    r = math.isqrt(n)
    return r * r == n, r


def sieve_primes(limit: int) -> List[int]:
    if limit < 2:
        return []

    is_prime = [True] * (limit + 1)
    is_prime[0:2] = [False, False]

    for p in range(2, math.isqrt(limit) + 1):
        if is_prime[p]:
            for multiple in range(p * p, limit + 1, p):
                is_prime[multiple] = False

    return [n for n in range(2, limit + 1) if is_prime[n]]


def is_probable_prime(n: int) -> bool:
    """Miller-Rabin primality test. Deterministic for the input sizes used here."""
    if n < 2:
        return False
    small_bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for p in small_bases:
        if n % p == 0:
            return n == p

    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2

    for a in small_bases:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def tonelli_shanks(n: int, p: int) -> Optional[int]:
    """
    Find t such that:

        t^2 ≡ n (mod p)

    where p is prime.

    Return one root t if it exists, otherwise return None.
    The other root is (-t) % p.
    """
    n %= p
    if n == 0:
        return 0
    if p == 2:
        return n
    if pow(n, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        t = pow(n, (p + 1) // 4, p)
        return t

    q = p - 1
    s = 0
    while q % 2 == 0:
        s += 1
        q //= 2

    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1

    # t ^ 2 = n * u (mod p), initially u = n and t = n^((q + 1) // 2).
    m = s
    c = pow(z, q, p)
    u = pow(n, q, p)
    t = pow(n, (q + 1) // 2, p)

    while u != 1:
        i = 1
        u2 = pow(u, 2, p)

        while u2 != 1:
            u2 = pow(u2, 2, p)
            i += 1
            if i == m:
                return None
        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = pow(b, 2, p)
        u = (u * c) % p
        t = (t * b) % p

    return t


def modular_square_roots_prime(N: int, p: int) -> List[int]:
    """Roots of x^2 = N (mod p)."""
    if p == 2:
        return [r for r in (0, 1) if (r * r - N) % 2 == 0]
    r = tonelli_shanks(N, p)
    if r is None:
        return []
    return sorted(set([r % p, (-r) % p]))


def estimate_parameters_task1(N: int, multiplier: float = 6.0, max_B: int = 50_000) -> QSParameters:
    """
    Estimate QS parameters using L(N) = exp(sqrt(ln N * ln ln N)).

    theoretical_B is a standard simple estimate exp(0.5 * sqrt(ln N ln ln N)).
    used_B is intentionally larger for this educational Python implementation so that
    90-100 bit examples usually finish within reasonable time.
    """
    if N <= 2:
        raise ValueError("N must be greater than 2")
    ln_n = math.log(N)
    ln_ln_n = math.log(ln_n)
    L_n = math.exp(math.sqrt(ln_n * ln_ln_n))
    theoretical_B = int(math.exp((1 / math.sqrt(2)) * math.sqrt(ln_n * ln_ln_n)))
    used_B = max(100, min(max_B, int(theoretical_B * multiplier)))
    return QSParameters(
        n_bits=N.bit_length(),
        ln_n=ln_n,
        ln_ln_n=ln_ln_n,
        L_n=L_n,
        theoretical_B=theoretical_B,
        used_B=used_B,
    )


def build_factor_base(N: int, B: int) -> Tuple[Optional[int], List[int], Dict[int, List[int]]]:
    """
    Build factor base: primes p <= B for which N is a quadratic residue mod p.

    Return (small_factor, factor_base, roots_by_prime). If p divides N, small_factor is returned.
    """
    factor_base: List[int] = []
    roots_by_prime: Dict[int, List[int]] = {}

    for p in sieve_primes(B):
        if N % p == 0:
            return p, [], {}
        roots = modular_square_roots_prime(N, p)
        if roots:
            factor_base.append(p)
            roots_by_prime[p] = roots

    return None, factor_base, roots_by_prime


def factor_over_base(f: int, factor_base: List[int]) -> Tuple[bool, List[int]]:
    """Exact integer verification that f is B-smooth over factor_base."""
    f_value = f
    exponents: List[int] = []

    for p in factor_base:
        count = 0
        while f_value % p == 0:
            f_value //= p
            count += 1
        exponents.append(count)
        if f_value == 1:
            exponents.extend([0] * (len(factor_base) - len(exponents)))
            return True, exponents

    return f_value == 1, exponents


def exponents_to_row_bits(exponents: List[int]) -> int:
    """Convert exponent vector to parity bit row over GF(2)."""
    bits = 0
    for i, exponent in enumerate(exponents):
        if exponent & 1:
            bits |= 1 << i
    return bits


def find_b_smooth_relations_task2(
    N: int,
    factor_base: List[int],
    roots_by_prime: Dict[int, List[int]],
    offset: int,
    length: int,
    threshold: Optional[float] = None,
) -> Tuple[List[Relation], int]:
    """
    Sieve f(t) = t^2 - N for x = ceil(sqrt(N)) + offset ... + offset + length - 1.

    For speed, the sieve subtracts log(p) at roots modulo p. Exact integer division
    is used afterwards, so false candidates do not affect correctness.
    """
    if np is None:  # pragma: no cover
        raise RuntimeError(f"NumPy is required for this implementation: {_NUMPY_IMPORT_ERROR}")

    a = ceil_sqrt(N)
    base = a + offset

    index = np.arange(length, dtype=np.float64)
    q0 = float(base * base - N) 
    #  (base + i) ** 2 - N = (base ** 2 - N) + 2 * base * i + i ** 2
    logs = np.log(q0 + 2.0 * float(base) * index + index * index)

    for p in factor_base:
        log_p = math.log(p)
        for root in roots_by_prime[p]:
            start = (root - base) % p
            logs[start::p] -= log_p

    if threshold is None:
        threshold = 2.0 * math.log(factor_base[-1])

    zero_indices = np.where(logs <= threshold)[0]
    relations: List[Relation] = []

    for i in zero_indices.tolist():
        x = base + i
        f_value = x * x - N
        is_smooth, exponents = factor_over_base(f_value, factor_base)
        if is_smooth:
            relations.append(
                Relation(
                    t=x,
                    f_value=f_value,
                    exponents=exponents,
                    row_bits=exponents_to_row_bits(exponents),
                )
            )

    return relations, int(zero_indices.size)


def find_dependencies_task3(row_bits: Iterable[int]) -> Iterable[int]:
    """
    Incremental Gaussian elimination over GF(2).

    Input rows are bitsets. A yielded integer is a bitset selecting relation rows whose XOR is zero.
    """
    pivots: Dict[int, Tuple[int, int]] = {}

    for i, row in enumerate(row_bits):
        result = row
        combination = 1 << i

        while result:
            pivot_col = result.bit_length() - 1
            if pivot_col in pivots:
                result ^= pivots[pivot_col][0]
                combination ^= pivots[pivot_col][1]
            else:
                pivots[pivot_col] = (result, combination)
                break

        if result == 0:
            yield combination


def try_dependency_task4(
    dependency: int,
    relations: List[Relation],
    factor_base: List[int],
    N: int,
) -> Optional[Tuple[int, int, int, int, int, int]]:
    """
    Given a dependency, build X^2 = Y^2 (mod N) and try gcd.

    Return (factor, X, Y, used_relations, gcd_minus, gcd_plus) or None.
    """
    X = 1
    exponent_sums = [0] * len(factor_base)
    used_count = 0

    index = 0
    dep = dependency
    while dep:
        if dep & 1:
            used_count += 1
            rel = relations[index]
            X = (X * rel.t) % N
            for j, exponent in enumerate(rel.exponents):
                exponent_sums[j] += exponent
        dep >>= 1
        index += 1

    Y = 1
    for p, exponent_sum in zip(factor_base, exponent_sums):
        if exponent_sum & 1:
            return None
        Y = (Y * pow(p, exponent_sum // 2, N)) % N

    gcd_minus = math.gcd((X - Y) % N, N)
    if 1 < gcd_minus < N:
        return gcd_minus, X, Y, used_count, gcd_minus, math.gcd((X + Y) % N, N)

    gcd_plus = math.gcd((X + Y) % N, N)
    if 1 < gcd_plus < N:
        return gcd_plus, X, Y, used_count, gcd_minus, gcd_plus

    return None


def quadratic_sieve(
    N: int,
    *,
    B: Optional[int] = None,
    max_rounds: int = 8,
    verbose: bool = True,
    save_json: Optional[str | Path] = None,
) -> int:
    """Return a non-trivial factor of composite N using Quadratic Sieve."""
    if N <= 1:
        raise ValueError("N must be greater than 1")
    if np is None:  # pragma: no cover
        raise RuntimeError(f"NumPy is required: {_NUMPY_IMPORT_ERROR}")

    follow = QSTrace(N=N)

    square, root = is_square(N)
    if square:
        return root
    if N % 2 == 0:
        return 2
    if is_probable_prime(N):
        return N

    params = estimate_parameters_task1(N)
    if B is not None:
        params.used_B = B
    follow.parameters = params

    current_B = params.used_B
    for attempt in range(4):
        small_factor, factor_base, roots_by_prime = build_factor_base(N, current_B)
        if small_factor is not None:
            follow.result = {
                "factor": small_factor,
                "other_factor": N // small_factor,
                "check_product": small_factor * (N // small_factor),
            }
            if save_json is not None:
                follow.save_json(save_json)
            return small_factor
        if not factor_base:
            raise RuntimeError("Empty factor base; choose a larger B")

        required_relations = len(factor_base) + 25
        sieve_length = max(20_000, current_B * 35)

        params.used_B = current_B
        params.factor_base_size = len(factor_base)
        params.required_relations = required_relations
        params.sieve_length = sieve_length
        follow.parameters = params
        follow.factor_base = factor_base

        if verbose:
            print(
                f"[parameters] bits={N.bit_length()} "
                f"L(N)≈{params.L_n:.3e} theoretical_B={params.theoretical_B} "
                f"used_B={current_B} factor_base_size={len(factor_base)} "
                f"required_relations≈{required_relations} sieve_length={sieve_length}",
                flush=True,
            )
            print(
                f"[factor base] first primes={factor_base[:15]}"
                f"{' ...' if len(factor_base) > 15 else ''}",
                flush=True,
            )

        relations: List[Relation] = []
        seen_x = set()
        offset = 0

        for round_no in range(1, max_rounds + 1):
            start_time = time.time()
            batch, candidates = find_b_smooth_relations_task2(N, factor_base, roots_by_prime, offset, sieve_length)
            new_count = 0
            for relation in batch:
                if relation.t not in seen_x:
                    seen_x.add(relation.t)
                    relations.append(relation)
                    new_count += 1

            round_info = {
                "round": round_no,
                "offset_start": offset,
                "offset_end": offset + sieve_length,
                "candidate_count": candidates,
                "new_smooth_relations": new_count,
                "total_smooth_relations": len(relations),
                "seconds": round(time.time() - start_time, 4),
            }
            follow.sieve_rounds.append(round_info)

            if verbose:
                print(
                    f"[sieve] round={round_no} interval=[{offset}, {offset + sieve_length}) "
                    f"candidates={candidates} new_smooth={new_count} "
                    f"total_smooth={len(relations)} time={round_info['seconds']}s",
                    flush=True,
                )

            if len(relations) >= required_relations:
                row_bits = [relation.row_bits for relation in relations]
                follow.matrix_rows_binary = [bin(row) for row in row_bits]
                follow.relations = [relation.to_light_dict() for relation in relations]

                if verbose:
                    print(
                        f"[linear algebra] solving GF(2) matrix with "
                        f"{len(relations)} rows and {len(factor_base)} columns",
                        flush=True,
                    )

                tried = 0
                for dependency in find_dependencies_task3(row_bits):
                    tried += 1
                    result = try_dependency_task4(dependency, relations, factor_base, N)
                    if result is not None:
                        factor, X, Y, used_count, gcd_minus, gcd_plus = result
                        follow.dependencies_tried += tried
                        follow.congruence = {
                            "X": X,
                            "Y": Y,
                            "used_relations": used_count,
                            "gcd_X_minus_Y_N": gcd_minus,
                            "gcd_X_plus_Y_N": gcd_plus,
                        }
                        follow.result = {
                            "factor": factor,
                            "other_factor": N // factor,
                            "check_product": factor * (N // factor),
                        }
                        if save_json is not None:
                            follow.save_json(save_json)

                        if verbose:
                            print(f"[linear algebra] dependencies_tried={tried}", flush=True)
                            print(f"[congruence] X={X}", flush=True)
                            print(f"[congruence] Y={Y}", flush=True)
                            print(f"[congruence] used_relations={used_count}", flush=True)
                            print(f"[gcd] gcd(X-Y, N)={gcd_minus}", flush=True)
                            print(f"[gcd] gcd(X+Y, N)={gcd_plus}", flush=True)
                        return factor

                if verbose:
                    print("[linear algebra] only trivial dependencies so far; collecting more relations", flush=True)
                required_relations += 25
                params.required_relations = required_relations

            offset += sieve_length

        current_B = int(current_B * 1.5) + 100
        if verbose:
            print(f"[retry] increasing B to {current_B}", flush=True)

    if save_json is not None:
        follow.save_json(save_json)
    raise RuntimeError("Quadratic Sieve failed; try larger --B or larger --max-rounds")


def factor_integer(
    N: int,
    *,
    verbose: bool = True,
    B: Optional[int] = None,
    save_json: Optional[str | Path] = None,
) -> List[int]:
    """Fully factor N recursively. The main non-trivial splitting method is QS."""
    if N == 1:
        return []
    if is_probable_prime(N):
        return [N]

    factor = quadratic_sieve(N, B=B, verbose=verbose, save_json=save_json)
    if factor == N:
        return [N]

    return sorted(
        factor_integer(factor, verbose=False, B=B)
        + factor_integer(N // factor, verbose=False, B=B)
    )
