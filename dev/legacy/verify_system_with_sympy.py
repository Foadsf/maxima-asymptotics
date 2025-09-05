#!/usr/bin/env python3
"""
Verification script for system ODE test cases using SymPy.
"""
import sympy as sp


def verify_system_case(name, ode_list, funcs, ics, order):
    """Verifies a single system ODE case."""
    print(f"\n{'='*60}")
    print(f"Case: {name}")
    print(f"{'='*60}")

    t = funcs[0].args[0]

    try:
        # Use dsolve for systems of ODEs
        exact_solutions = sp.dsolve(ode_list, funcs, ics=ics)
        print("Exact solutions:", exact_solutions)

        series_polys = []
        for sol in exact_solutions:
            series_exp = sp.series(sol.rhs, t, 0, order + 1)
            series_polys.append(series_exp.removeO())

        print("\nExpected Taylor Series Polynomials:")
        for func, poly in zip(funcs, series_polys):
            print(f"  {func}: {poly}")

        return series_polys

    except Exception as e:
        print(f"Could not solve the system exactly: {e}")
        # Manual calculation for simple systems
        print("Attempting manual derivative calculation...")
        try:
            # This is complex to generalize. For now, we state the expected result.
            print(
                "Manual verification is complex; rely on known series for this script."
            )
        except:
            pass
        return None


def main():
    t = sp.symbols("t")
    f = sp.Function("f")(t)
    g = sp.Function("g")(t)

    # Test cases from the Maxima suite
    cases = [
        (
            "sin/cos Oscillator",
            [sp.Eq(f.diff(t), g), sp.Eq(g.diff(t), -f)],
            [f, g],
            {f.subs(t, 0): 0, g.subs(t, 0): 1},
            5,
        ),
        (
            "Exponential System",
            [sp.Eq(f.diff(t), f + g), sp.Eq(g.diff(t), f + g)],
            [f, g],
            {f.subs(t, 0): 1, g.subs(t, 0): 0},
            5,
        ),
        (
            "Mixed Polynomial System",
            [sp.Eq(f.diff(t), t), sp.Eq(g.diff(t), f + g)],
            [f, g],
            {f.subs(t, 0): 0, g.subs(t, 0): 1},
            4,
        ),
    ]

    for name, odes, funcs, ics, order in cases:
        verify_system_case(name, odes, funcs, ics, order)

    print(f"\n{'='*60}")
    print("Compare these results with your Maxima test suite expectations.")


if __name__ == "__main__":
    main()
