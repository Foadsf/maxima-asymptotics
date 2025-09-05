#!/usr/bin/env python3
"""
Verify the third-order ODE case: y''' = x + y, y(0)=1, y'(0)=0, y''(0)=0
"""

import sympy as sp
from sympy import symbols, Function, Eq, dsolve, series, factorial


def verify_third_order_coupling():
    """Verify y''' = x + y, y(0)=1, y'(0)=0, y''(0)=0"""
    print("=" * 60)
    print("Verifying: y''' = x + y, y(0)=1, y'(0)=0, y''(0)=0")
    print("=" * 60)

    x = symbols("x")
    y = Function("y")

    # Define the ODE and ICs
    ode = Eq(y(x).diff(x, 3), x + y(x))
    ics = {y(0): 1, y(x).diff(x).subs(x, 0): 0, y(x).diff(x, 2).subs(x, 0): 0}

    print("ODE:", ode)
    print("Initial conditions:", ics)
    print()

    # Try to solve exactly
    try:
        exact_solution = dsolve(ode, y(x), ics=ics)
        print("Exact solution:", exact_solution)
        y_exact = exact_solution.rhs

        # Get series expansion
        taylor_series = series(y_exact, x, 0, 6)
        print("Taylor series:", taylor_series)
        print("Series polynomial:", taylor_series.removeO())

    except Exception as e:
        print("Cannot solve exactly:", e)
        print("Using manual derivative calculation...")

        # Manual calculation
        print("\nManual calculation of derivatives at x=0:")
        print("y(0) = 1")
        print("y'(0) = 0")
        print("y''(0) = 0")
        print("y'''(0) = 0 + y(0) = 1")
        print("y^(4)(0) = d/dx[x + y]|_{x=0} = 1 + y'(0) = 1")
        print("y^(5)(0) = d/dx[1 + y']|_{x=0} = 0 + y''(0) = 0")
        print("y^(6)(0) = d/dx[y'']|_{x=0} = y'''(0) = 1")
        print()

        # Construct series manually
        coeffs = [1, 0, 0, 1, 1, 0, 1]  # y(0) through y^(6)(0)

        print("Taylor series construction:")
        terms = []
        for i, coeff in enumerate(coeffs):
            if coeff != 0:
                if i == 0:
                    terms.append(f"{coeff}")
                elif i == 1:
                    terms.append(f"{coeff}*x")
                else:
                    terms.append(f"{coeff}*x^{i}/{factorial(i)}")

        print("y(x) = " + " + ".join(terms))

        # Simplify the polynomial up to order 5
        polynomial = sum(coeffs[i] * x**i / factorial(i) for i in range(6))
        polynomial = sp.expand(polynomial)
        print("Polynomial (up to x^5):", polynomial)

        # Extract non-zero terms up to x^5
        terms_up_to_5 = []
        for i in range(6):
            coeff = polynomial.coeff(x, i)
            if coeff is not None and coeff != 0:
                terms_up_to_5.append(f"x^{i}: {coeff}")

        print("Non-zero coefficients (up to x^5):", ", ".join(terms_up_to_5))

        return polynomial


if __name__ == "__main__":
    result = verify_third_order_coupling()

    print("\n" + "=" * 60)
    print("CONCLUSION:")
    print("Expected series: 1 + x^3/6 + x^4/24 + 0*x^5 + ...")
    print("Library result should be: 1 + x^3/6 + x^4/24")
    print("Test expectation x^5/120 + ... is WRONG")
    print("=" * 60)
