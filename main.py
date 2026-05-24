from __future__ import annotations

import argparse
import math
import time

from qs_course.quadratic_sieve import factor_integer, quadratic_sieve


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Mini-course: integer factorization by Quadratic Sieve"
    )
    parser.add_argument("N", type=int, help="Composite integer to factor")
    parser.add_argument(
        "--B",
        type=int,
        default=None,
        help="Optional factor-base bound. If omitted, the program estimates it from N.",
    )
    parser.add_argument(
        "--max-rounds",
        type=int,
        default=8,
        help="Maximum sieve intervals before increasing B.",
    )
    parser.add_argument(
        "--one-factor",
        action="store_true",
        help="Return only one non-trivial factor instead of recursive full factorization.",
    )
    parser.add_argument(
        "--save-json",
        default=None,
        help="Optional path for saving intermediate QS data as JSON.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    N = args.N
    start = time.time()

    print(f"[input] N={N}")
    print(f"[input] bit_length={N.bit_length()}")

    if args.one_factor:
        factor = quadratic_sieve(
            N,
            B=args.B,
            max_rounds=args.max_rounds,
            verbose=True,
            save_json=args.save_json,
        )
        other = N // factor
        print(f"[result] non_trivial_factor={factor}")
        print(f"[result] other_factor={other}")
        print(f"[check] {factor} * {other} = {factor * other}")
    else:
        factors = factor_integer(N, verbose=True, B=args.B, save_json=args.save_json)
        print(f"[result] factors={factors}")
        print(f"[result] {' * '.join(map(str, factors))} = {math.prod(factors)}")
        print(f"[check] product_equals_N={math.prod(factors) == N}")

    print(f"[time] {time.time() - start:.3f}s")


if __name__ == "__main__":
    main()
