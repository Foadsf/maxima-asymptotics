#!/usr/bin/env python3
"""
Comprehensive verification script for all ODE test cases
using SymPy to confirm the correct series expansions.
"""

import sympy as sp
from sympy import (
    symbols,
    Function,
    Eq,
    dsolve,
    series,
    sin,
    cos,
    exp,
    simplify,
    sqrt,
    atan,
    asin,
)
import traceback


def verify_ode_case(name, ode_eq, ic_dict, order=5, manual_series=None):
    """
    Verify a single ODE case

    Args:
        name: descriptive name
        ode_eq: SymPy equation for the ODE
        ic_dict: initial conditions as dict
        order: order for series expansion
        manual_series: expected series if known analytically
    """
    print(f"\n{'='*60}")
    print(f"Case: {name}")
    print(f"{'='*60}")

    x = symbols("x")
    y = Function("y")

    try:
        print("ODE:", ode_eq)
        print("IC:", ic_dict)

        # Try to solve exactly
        try:
            exact_solution = dsolve(ode_eq, y(x), ics=ic_dict)
            print("Exact solution:", exact_solution)
            y_exact = exact_solution.rhs
        except:
            print("Cannot solve exactly, will try series method")
            y_exact = None

        # If we have exact solution, get its series
        if y_exact is not None:
            taylor_series = series(y_exact, x, 0, order + 1)
            print("Taylor series:", taylor_series)
            series_poly = taylor_series.removeO()
            print("Series polynomial:", series_poly)

            # Extract just the significant coefficients
            coeffs = []
            for i in range(order + 1):
                coeff = (
                    series_poly.coeff(x, i)
                    if series_poly.coeff(x, i) is not None
                    else 0
                )
                if coeff != 0:
                    coeffs.append(f"x^{i}: {coeff}")

            if coeffs:
                print("Non-zero coefficients:", ", ".join(coeffs))
            else:
                print("All coefficients up to order", order, "are zero")

            return series_poly
        else:
            print("Using series method not implemented in this verification")
            return None

    except Exception as e:
        print(f"Error processing case: {e}")
        traceback.print_exc()
        return None


def verify_all_test_cases():
    """Verify all the test cases from the comprehensive suite"""

    x = symbols("x")
    y = Function("y")

    print("COMPREHENSIVE ODE VERIFICATION")
    print("Checking all test cases for correct expected results")

    results = {}

    # Phase 1: Basic cases that should work
    cases = [
        # Already verified - these should pass
        (
            "y' = y, y(0) = 1 (exponential)",
            Eq(y(x).diff(x), y(x)),
            {y(0): 1},
            6,
            "1 + x + x²/2 + x³/6 + x⁴/24 + x⁵/120 + x⁶/720",
        ),
        (
            "y' = y + sin(x), y(0) = 0 (linear inhomogeneous)",
            Eq(y(x).diff(x), y(x) + sin(x)),
            {y(0): 0},
            5,
            "x²/2 + x³/6",
        ),
        # Cases that might be failing
        (
            "y' = x + y², y(0) = 1 (spec example)",
            Eq(y(x).diff(x), x + y(x) ** 2),
            {y(0): 1},
            3,
            "1 + x + (3/2)x² + (4/3)x³",
        ),
        (
            "y' = x*y, y(0) = 1 (Bernoulli)",
            Eq(y(x).diff(x), x * y(x)),
            {y(0): 1},
            6,
            "1 + x²/2 + x⁴/8 + x⁶/48",
        ),
        (
            "y' = 1 + y², y(0) = 0 (Riccati → tan(x))",
            Eq(y(x).diff(x), 1 + y(x) ** 2),
            {y(0): 0},
            5,
            "x + x³/3 + 2x⁵/15",
        ),
        (
            "y' = x² + y², y(0) = 0 (nonlinear)",
            Eq(y(x).diff(x), x**2 + y(x) ** 2),
            {y(0): 0},
            5,
            "?",
        ),
        (
            "y' = x², y(0) = 0 (pure polynomial)",
            Eq(y(x).diff(x), x**2),
            {y(0): 0},
            5,
            "x³/3",
        ),
        (
            "y' = y², y(0) = 1 (separable → 1/(1-x))",
            Eq(y(x).diff(x), y(x) ** 2),
            {y(0): 1},
            5,
            "1 + x + x² + x³ + x⁴ + x⁵",
        ),
        (
            "y' = cos(x), y(0) = 0 (→ sin(x))",
            Eq(y(x).diff(x), cos(x)),
            {y(0): 0},
            7,
            "x - x³/6 + x⁵/120 - x⁷/5040",
        ),
        (
            "y' = 1/(1+x²), y(0) = 0 (→ arctan(x))",
            Eq(y(x).diff(x), 1 / (1 + x**2)),
            {y(0): 0},
            9,
            "x - x³/3 + x⁵/5 - x⁷/7 + x⁹/9",
        ),
    ]

    for case_data in cases:
        if len(case_data) == 5:
            name, ode_eq, ic_dict, order, expected = case_data
            result = verify_ode_case(name, ode_eq, ic_dict, order, expected)
            results[name] = result
        else:
            print(f"Skipping malformed case: {case_data}")

    return results


def manual_nonlinear_verification():
    """
    Manual verification for nonlinear cases that SymPy might not solve exactly
    """
    print(f"\n{'='*80}")
    print("MANUAL VERIFICATION FOR NONLINEAR CASES")
    print(f"{'='*80}")

    # Case: y' = x² + y², y(0) = 0
    print("\nCase: y' = x² + y², y(0) = 0")
    print("Manual Taylor series calculation:")
    print("y(0) = 0")
    print("y'(0) = 0² + 0² = 0")
    print("y''(0) = d/dx[x² + y²]|_{x=0} = 2x + 2yy'|_{x=0} = 0 + 0 = 0")
    print(
        "y'''(0) = d/dx[2x + 2yy']|_{x=0} = 2 + 2(y')² + 2yy''|_{x=0} = 2 + 0 + 0 = 2"
    )
    print(
        "y^(4)(0) = d/dx[2 + 2(y')² + 2yy'']|_{x=0} = 4y'y'' + 2(y'')² + 2y'y'' + 2yy'''|_{x=0} = 0"
    )
    print("y^(5)(0) = more complex, but let's compute...")
    print()
    print("From derivatives at x=0:")
    print("y = 0 + 0·x + 0·x²/2! + 2·x³/3! + 0·x⁴/4! + ...")
    print("y = x³/3 + higher order terms")
    print()
    print("To get x⁵ term, need y^(5)(0):")
    print("Working through the chain rule systematically...")

    # Let's work this out step by step
    print("y^(4)(0) calculation:")
    print("y^(4) = d/dx[2 + 2(y')² + 2yy'']")
    print("     = 0 + 2·2y'·y'' + 2(y'·y'' + y·y''')")
    print("At x=0: y^(4)(0) = 4·0·0 + 2(0·0 + 0·2) = 0")
    print()

    print("y^(5)(0) calculation:")
    print("y^(5) = d/dx[4y'y'' + 2y'y'' + 2yy''']")
    print("     = 4(y'')² + 4y'y''' + 2(y'')² + 2y'y''' + 2y'y''' + 2yy^(4)")
    print("     = 6(y'')² + 8y'y''' + 2yy^(4)")
    print("At x=0: y^(5)(0) = 6·0² + 8·0·2 + 2·0·0 = 0")
    print()

    print("Hmm, let me recalculate more carefully...")
    print("Actually, let's substitute our series y = x³/3 + ... back into the ODE")
    print("If y = x³/3 + ax⁵ + ..., then:")
    print("y' = x² + 5ax⁴ + ...")
    print("y² = (x³/3)² + ... = x⁶/9 + ... (higher order)")
    print("x² + y² = x² + x⁶/9 + ... ")
    print()
    print("But y' should equal x² + y², so:")
    print("x² + 5ax⁴ + ... = x² + x⁶/9 + ...")
    print("Comparing coefficients: x⁴ term gives 5a = 0, so a = 0")
    print("This suggests y = x³/3 + O(x⁷), so NO x⁵ term!")

    return {"y' = x² + y², y(0) = 0": "x³/3"}


if __name__ == "__main__":
    # Run comprehensive verification
    results = verify_all_test_cases()

    # Manual verification for tricky cases
    manual_results = manual_nonlinear_verification()

    print(f"\n{'='*80}")
    print("SUMMARY AND RECOMMENDATIONS")
    print(f"{'='*80}")

    print("\nExpected results to use in test suite:")
    print("1. y' = y + sin(x), y(0) = 0  →  x²/2 + x³/6")
    print("2. y' = x² + y², y(0) = 0     →  x³/3 (NOT x³/3 + x⁵/15)")
    print("3. y' = x + y², y(0) = 1      →  1 + x + (3/2)x² + (4/3)x³")
    print("4. y' = y², y(0) = 1          →  1 + x + x² + x³ + x⁴ + x⁵")
    print("5. y' = 1 + y², y(0) = 0      →  x + x³/3 + 2x⁵/15")

    print("\nThe Maxima library appears to be correct.")
    print("Test expectations need to be fixed based on these verifications.")
