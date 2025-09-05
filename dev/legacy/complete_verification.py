#!/usr/bin/env python3
"""
Complete verification of all test cases in comprehensive_test_suite.mac
Using SymPy to compute exact results, no manual calculations or guessing.
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
    tan,
)
import traceback


def solve_and_expand(name, ode_eq, ic_dict, order=5, method="exact"):
    """
    Solve ODE and get Taylor expansion using SymPy only

    Returns:
        (success, series_polynomial, error_msg)
    """
    x = symbols("x")
    y = Function("y")

    try:
        if method == "exact":
            # Try exact solution first
            try:
                exact_solution = dsolve(ode_eq, y(x), ics=ic_dict)
                y_exact = exact_solution.rhs
                taylor_series = series(y_exact, x, 0, order + 1)
                return True, taylor_series.removeO(), None
            except Exception as e:
                return False, None, f"Could not solve exactly: {e}"

        elif method == "series":
            # For nonlinear ODEs, use SymPy's series method if available
            # This is a placeholder - SymPy doesn't have a direct series ODE solver
            return False, None, "Series method not implemented"

    except Exception as e:
        return False, None, f"Error: {e}"


def compute_derivatives_manually(ode_rhs, ic_values, order=5):
    """
    Compute derivatives manually by repeated differentiation
    ode_rhs: right-hand side of y' = f(x,y) or y''' = f(x,y,y',y'')
    ic_values: dict with y(0), y'(0), etc.
    """
    x, y_var = symbols("x y")

    # Replace y with a symbol for differentiation
    f = ode_rhs.subs(Function("y")(x), y_var)

    derivatives = {}
    derivatives[0] = ic_values.get("y", 0)

    # For first-order ODE y' = f(x,y)
    if 1 not in ic_values:  # y'(0) not given, compute from ODE
        derivatives[1] = f.subs([(x, 0), (y_var, derivatives[0])])
    else:
        derivatives[1] = ic_values[1]

    # Compute higher derivatives by chain rule
    current_f = f
    for k in range(2, order + 1):
        if k in ic_values:
            derivatives[k] = ic_values[k]
        else:
            # Differentiate current_f
            df_dx = sp.diff(current_f, x)
            df_dy = sp.diff(current_f, y_var)

            # Apply chain rule: d/dx[f(x,y)] = df/dx + df/dy * dy/dx
            next_f = (
                df_dx + df_dy * derivatives[1]
            )  # This is approximate for higher order

            # Evaluate at x=0 with known derivative values
            subs_list = [(x, 0), (y_var, derivatives[0])]
            derivatives[k] = next_f.subs(subs_list)
            current_f = next_f

    # Build polynomial
    polynomial = sum(derivatives[k] * x**k / sp.factorial(k) for k in range(order + 1))
    return sp.expand(polynomial), derivatives


def verify_all_cases():
    """Systematically verify all test cases"""
    x = symbols("x")
    y = Function("y")

    test_cases = [
        # Phase 1: Basic cases
        ("y' = y, y(0) = 1", Eq(y(x).diff(x), y(x)), {y(0): 1}, 6),
        ("y' = 2*y, y(0) = 3", Eq(y(x).diff(x), 2 * y(x)), {y(0): 3}, 4),
        ("y' = x^2, y(0) = 0", Eq(y(x).diff(x), x**2), {y(0): 0}, 5),
        ("y' = x + x^2, y(0) = 1", Eq(y(x).diff(x), x + x**2), {y(0): 1}, 4),
        ("y' = y^2, y(0) = 1", Eq(y(x).diff(x), y(x) ** 2), {y(0): 1}, 5),
        # Phase 1: Nonlinear cases
        ("y' = x + y^2, y(0) = 1", Eq(y(x).diff(x), x + y(x) ** 2), {y(0): 1}, 3),
        ("y' = x*y, y(0) = 1", Eq(y(x).diff(x), x * y(x)), {y(0): 1}, 6),
        ("y' = y + sin(x), y(0) = 0", Eq(y(x).diff(x), y(x) + sin(x)), {y(0): 0}, 5),
        ("y' = 1 + y^2, y(0) = 0", Eq(y(x).diff(x), 1 + y(x) ** 2), {y(0): 0}, 5),
        ("y' = x^2 + y^2, y(0) = 0", Eq(y(x).diff(x), x**2 + y(x) ** 2), {y(0): 0}, 5),
        # Non-zero expansion points
        ("y' = x^2, y(1) = 2", Eq(y(x).diff(x), x**2), {y(1): 2}, 4),
        ("y' = y, y(2) = 1", Eq(y(x).diff(x), y(x)), {y(2): 1}, 4),
        ("y' = 2*x, y(-1) = 0", Eq(y(x).diff(x), 2 * x), {y(-1): 0}, 3),
        # Edge cases
        ("y' = 5, y(0) = 1", Eq(y(x).diff(x), 5), {y(0): 1}, 3),
        # Special functions
        ("y' = cos(x), y(0) = 0", Eq(y(x).diff(x), cos(x)), {y(0): 0}, 7),
        ("y' = exp(x), y(0) = 1", Eq(y(x).diff(x), exp(x)), {y(0): 1}, 5),
        ("y' = 1/(1+x^2), y(0) = 0", Eq(y(x).diff(x), 1 / (1 + x**2)), {y(0): 0}, 9),
        # Phase 2: Second-order ODEs
        (
            "y'' = y, y(0)=1, y'(0)=0",
            Eq(y(x).diff(x, 2), y(x)),
            {y(0): 1, y(x).diff(x).subs(x, 0): 0},
            6,
        ),
        (
            "y'' = -y, y(0)=1, y'(0)=0",
            Eq(y(x).diff(x, 2), -y(x)),
            {y(0): 1, y(x).diff(x).subs(x, 0): 0},
            6,
        ),
        (
            "y'' = -y, y(0)=0, y'(0)=1",
            Eq(y(x).diff(x, 2), -y(x)),
            {y(0): 0, y(x).diff(x).subs(x, 0): 1},
            5,
        ),
        (
            "y'' = -4*y, y(0)=1, y'(0)=0",
            Eq(y(x).diff(x, 2), -4 * y(x)),
            {y(0): 1, y(x).diff(x).subs(x, 0): 0},
            4,
        ),
        # Phase 2: Higher-order ODEs
        (
            "y''' = 0, y(0)=1, y'(0)=2, y''(0)=3",
            Eq(y(x).diff(x, 3), 0),
            {y(0): 1, y(x).diff(x).subs(x, 0): 2, y(x).diff(x, 2).subs(x, 0): 3},
            3,
        ),
        (
            "y^(4) = y, y(0)=1, y'(0)=0, y''(0)=0, y'''(0)=0",
            Eq(y(x).diff(x, 4), y(x)),
            {
                y(0): 1,
                y(x).diff(x).subs(x, 0): 0,
                y(x).diff(x, 2).subs(x, 0): 0,
                y(x).diff(x, 3).subs(x, 0): 0,
            },
            8,
        ),
        (
            "y''' = x + y, y(0)=1, y'(0)=0, y''(0)=0",
            Eq(y(x).diff(x, 3), x + y(x)),
            {y(0): 1, y(x).diff(x).subs(x, 0): 0, y(x).diff(x, 2).subs(x, 0): 0},
            5,
        ),
    ]

    results = {}

    print("COMPLETE VERIFICATION OF ALL TEST CASES")
    print("=" * 60)

    for case_name, ode, ics, order in test_cases:
        print(f"\nCase: {case_name}")
        print("-" * 40)

        # Try exact solution
        success, result, error = solve_and_expand(case_name, ode, ics, order)

        if success:
            print(f"SUCCESS: {result}")
            results[case_name] = result
        else:
            print(f"FAILED: {error}")

            # For nonlinear cases that can't be solved exactly,
            # try numerical coefficient computation via derivatives
            if "y^2" in case_name or "x + y^2" in case_name or "x^2 + y^2" in case_name:
                print("Attempting manual derivative computation...")
                try:
                    # Extract RHS and convert to computable form
                    rhs = ode.rhs
                    # This is complex for higher-order, skip for now
                    print("Manual computation not implemented for this case")
                except:
                    print("Could not compute manually")

            results[case_name] = None

    return results


def format_for_maxima(expr):
    """Convert SymPy expression to Maxima-friendly format"""
    if expr is None:
        return "UNKNOWN"

    # Convert SymPy to string and make Maxima-compatible
    maxima_str = str(expr)
    maxima_str = maxima_str.replace("**", "^")
    return maxima_str


def generate_test_expectations(results):
    """Generate the corrected expectations for the test suite"""
    print("\n" + "=" * 80)
    print("CORRECTED EXPECTATIONS FOR TEST SUITE")
    print("=" * 80)

    for case_name, result in results.items():
        if result is not None:
            maxima_format = format_for_maxima(result)
            print(f"{case_name:40} -> {maxima_format}")
        else:
            print(f"{case_name:40} -> COULD NOT VERIFY")

    print("\n" + "=" * 80)
    print("Use these results to update comprehensive_test_suite.mac")
    print("=" * 80)


if __name__ == "__main__":
    results = verify_all_cases()
    generate_test_expectations(results)
