#!/usr/bin/env python3
"""
sympy_crosscheck_all.py

Cross-checks all hardcoded assumptions (expected Taylor polynomials) from
Maxima's comprehensive_test_suite.mac using SymPy.

- Uses exact solutions (dsolve) when available.
- Falls back to a robust, generic power-series recursion for:
    * First-order ODEs:          y' = F(x, y)
    * m-th order ODEs:           y^(m) = G(x, y, y', ..., y^(m-1))
    * 2×2 first-order systems:   f' = F(x, f, g), g' = G(x, f, g)

Outputs a PASS/FAIL report and prints the SymPy-derived series alongside the
Maxima-expected series for any mismatches, plus suggested fixes.

Run:
    python3 sympy_crosscheck_all.py
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Union
import sympy as sp

# ----------------------------------------------------------------------------
# Utilities
# ----------------------------------------------------------------------------

x = sp.Symbol("x")
I = sp.I  # imaginary unit


def poly_to_maxima_str(expr: sp.Expr, var: sp.Symbol = x) -> str:
    """Pretty-print a polynomial/series for easy copy into Maxima (use '^')."""
    s = sp.simplify(sp.expand(expr))
    return str(s).replace("**", "^")


def series_from_exact(
    ode: Union[sp.Eq, Sequence[sp.Eq]], func, ics: Dict, order: int, x0=0
) -> Optional[sp.Expr]:
    """
    Try an exact solve (dsolve) and return the Taylor polynomial about x0 to 'order'.
    - ode: SymPy Eq or list of Eq for systems
    - func: y(x) for scalar, or [f(t), g(t)] for systems
    - ics: dictionary of initial conditions
    - order: truncation degree
    - x0: expansion point
    Returns:
        For scalar ODEs: sp.Expr polynomial in (x-x0)
        For systems: List[sp.Expr] polynomials in (x-x0)
    """
    try:
        sol = sp.dsolve(ode, func, ics=ics)
        if isinstance(sol, list):
            polynomials = []
            for s in sol:
                ser = sp.series(sp.expand(s.rhs), x, x0, order + 1).removeO()
                polynomials.append(sp.expand(ser))
            return polynomials  # type: ignore
        else:
            ser = sp.series(sp.expand(sol.rhs), x, x0, order + 1).removeO()
            return sp.expand(ser)
    except Exception:
        return None


# ----------------------------------------------------------------------------
# Generic power-series solvers
# ----------------------------------------------------------------------------


def series_first_order(F: sp.Expr, y0: sp.Expr, order: int, x0: sp.Expr = 0) -> sp.Expr:
    """
    Compute series for y' = F(x, y), y(x0) = y0 up to 'order' in (x-x0).
    Returns a polynomial in (x-x0).
    """
    t = sp.Symbol("t")
    # We'll expand around t = 0 where x = x0 + t
    # y(t) = sum_{k=0..order} a_k t^k
    a = [None] * (order + 1)  # coefficients
    a[0] = sp.sympify(y0)

    for n in range(order):  # we can determine a[n+1]
        Yn = sum((a[k] if a[k] is not None else 0) * t**k for k in range(n + 1))
        F_sub = sp.expand(
            F.subs({x: x0 + t, sp.Function("y")(x): Yn, sp.Symbol("y"): Yn})
        )
        F_ser = sp.series(F_sub, t, 0, n + 1).removeO()
        c_n = sp.expand(F_ser).coeff(t, n)  # coefficient of t^n in RHS
        a[n + 1] = sp.simplify(c_n / (n + 1))

    poly_t = sum((a[k] if a[k] is not None else 0) * t**k for k in range(order + 1))
    return sp.expand(poly_t.subs(t, x - x0))


def series_nth_order(
    m: int, G: sp.Expr, ics: Dict[int, sp.Expr], order: int, x0: sp.Expr = 0
) -> sp.Expr:
    """
    Compute series for y^(m) = G(x, y, y', ..., y^(m-1)) about x0 up to 'order'.
    ics: maps derivative order k (0..m-1) to y^(k)(x0).
    Returns polynomial in (x-x0).
    """
    t = sp.Symbol("t")
    # y(t) = sum_{k=0..order} a_k t^k
    a = [sp.S(0)] * (order + 1)
    for k, val in ics.items():
        a[k] = sp.sympify(val)

    # Helper to build truncated y and its derivatives using currently known a's
    def Y_upto(K: int):
        return sum(a[i] * t**i for i in range(min(K, order) + 1))

    for n in range(order - m + 1):
        # Build y, y', ..., y^(m-1) using known coefficients up to index n+m-1
        Y = Y_upto(n + m - 1)
        derivs = [sp.diff(Y, t, k) for k in range(m)]
        # Left side y^(m)
        # Right side G(x, y, y', ..., y^(m-1))
        subs_map = {x: x0 + t, sp.Function("y")(x): Y, sp.Symbol("y"): Y}
        for k in range(1, m):
            subs_map[sp.Function("y")(x).diff(x, k - 1)] = derivs[k - 1]

        G_sub = sp.expand(G.subs(subs_map))
        G_ser = sp.series(G_sub, t, 0, n + 1).removeO()
        c_n = sp.expand(G_ser).coeff(t, n)  # RHS coeff of t^n

        # Coefficient of t^n in y^(m) is a_{n+m} * (n+m)!/n!
        factor = sp.factorial(n + m) / sp.factorial(n)
        a[n + m] = sp.simplify(c_n / factor)

    poly_t = sum(a[k] * t**k for k in range(order + 1))
    return sp.expand(poly_t.subs(t, x - x0))


def series_system_2x2(
    F: sp.Expr, G: sp.Expr, f0: sp.Expr, g0: sp.Expr, order: int, x0: sp.Expr = 0
) -> Tuple[sp.Expr, sp.Expr]:
    """
    Compute series for the system:
        f' = F(x, f, g)
        g' = G(x, f, g)
    with f(x0)=f0, g(x0)=g0 up to 'order' in (x-x0).
    Returns (f_poly, g_poly).
    """
    t = sp.Symbol("t")
    af = [None] * (order + 1)
    ag = [None] * (order + 1)
    af[0] = sp.sympify(f0)
    ag[0] = sp.sympify(g0)

    for n in range(order):
        F_sub = sp.expand(
            F.subs(
                {
                    x: x0 + t,
                    sp.Function("f")(x): sum((af[k] or 0) * t**k for k in range(n + 1)),
                    sp.Function("g")(x): sum((ag[k] or 0) * t**k for k in range(n + 1)),
                }
            )
        )
        G_sub = sp.expand(
            G.subs(
                {
                    x: x0 + t,
                    sp.Function("f")(x): sum((af[k] or 0) * t**k for k in range(n + 1)),
                    sp.Function("g")(x): sum((ag[k] or 0) * t**k for k in range(n + 1)),
                }
            )
        )
        F_ser = sp.series(F_sub, t, 0, n + 1).removeO()
        G_ser = sp.series(G_sub, t, 0, n + 1).removeO()
        cF = sp.expand(F_ser).coeff(t, n)
        cG = sp.expand(G_ser).coeff(t, n)
        af[n + 1] = sp.simplify(cF / (n + 1))
        ag[n + 1] = sp.simplify(cG / (n + 1))

    poly_f = sum((af[k] or 0) * t**k for k in range(order + 1))
    poly_g = sum((ag[k] or 0) * t**k for k in range(order + 1))
    return sp.expand(poly_f.subs(t, x - x0)), sp.expand(poly_g.subs(t, x - x0))


# ----------------------------------------------------------------------------
# Test definitions mirroring comprehensive_test_suite.mac expectations
# ----------------------------------------------------------------------------


@dataclass
class ScalarCase:
    label: str
    ode_type: str  # 'first', 'nth'
    # For first-order: F(x, y)
    F: Optional[sp.Expr] = None
    # For nth-order: m, G(x, y, y', ..., y^(m-1))
    m: Optional[int] = None
    G: Optional[sp.Expr] = None
    # Initial conditions
    x0: sp.Expr = 0
    ics_first: Optional[sp.Expr] = None  # y(x0) = ...
    ics_nth: Optional[Dict[int, sp.Expr]] = None  # {k: y^(k)(x0)}
    order: int = 5
    expected: Optional[sp.Expr] = None  # polynomial in (x-x0)


@dataclass
class SystemCase:
    label: str
    F: sp.Expr
    G: sp.Expr
    x0: sp.Expr
    f0: sp.Expr
    g0: sp.Expr
    order: int
    expected_f: Optional[sp.Expr]
    expected_g: Optional[sp.Expr]


def build_cases() -> Tuple[List[ScalarCase], List[SystemCase]]:
    # Phase 1: Basic first-order
    P1_basic = [
        ScalarCase(
            "Linear homogeneous y'=y",
            "first",
            F=sp.Function("y")(x),
            x0=0,
            ics_first=1,
            order=6,
            expected=1
            + x
            + x**2 / sp.Integer(2)
            + x**3 / sp.Integer(6)
            + x**4 / sp.Integer(24)
            + x**5 / sp.Integer(120)
            + x**6 / sp.Integer(720),
        ),
        ScalarCase(
            "Linear y'=2y",
            "first",
            F=2 * sp.Function("y")(x),
            x0=0,
            ics_first=3,
            order=4,
            expected=3 + 6 * x + 6 * x**2 + 4 * x**3 + 2 * x**4,
        ),
        ScalarCase(
            "Polynomial source y'=x^2",
            "first",
            F=x**2,
            x0=0,
            ics_first=0,
            order=5,
            expected=x**3 / 3,
        ),
        ScalarCase(
            "Mixed polynomial y'=x+x^2",
            "first",
            F=x + x**2,
            x0=0,
            ics_first=1,
            order=4,
            expected=1 + x**2 / 2 + x**3 / 3,
        ),
        ScalarCase(
            "Separable y'=y^2",
            "first",
            F=sp.Function("y")(x) ** 2,
            x0=0,
            ics_first=1,
            order=5,
            expected=1 + x + x**2 + x**3 + x**4 + x**5,
        ),
    ]

    # Phase 1: Nonlinear first-order (only those with explicit expectations)
    P1_nonlinear = [
        ScalarCase(
            "Bernoulli y'=xy",
            "first",
            F=x * sp.Function("y")(x),
            x0=0,
            ics_first=1,
            order=6,
            expected=1 + x**2 / 2 + x**4 / 8 + x**6 / 48,
        ),
        ScalarCase(
            "Linear inhomo y'=y+sin(x)",
            "first",
            F=sp.Function("y")(x) + sp.sin(x),
            x0=0,
            ics_first=0,
            order=5,
            expected=x**2 / 2 + x**3 / 6,
        ),
        ScalarCase(
            "Riccati y'=1+y^2",
            "first",
            F=1 + sp.Function("y")(x) ** 2,
            x0=0,
            ics_first=0,
            order=5,
            expected=x + x**3 / 3 + 2 * x**5 / 15,
        ),
    ]

    # Phase 1: Non-zero expansion points
    t = x - 1
    P1_nonzero = [
        ScalarCase(
            "Nonzero point y'=x^2 at x=1",
            "first",
            F=x**2,
            x0=1,
            ics_first=2,
            order=4,
            expected=2 + t + t**2 + t**3 / sp.Integer(3),
        ),
        ScalarCase(
            "Exponential at x=2",
            "first",
            F=sp.Function("y")(x),
            x0=2,
            ics_first=1,
            order=4,
            expected=1
            + (x - 2)
            + (x - 2) ** 2 / 2
            + (x - 2) ** 3 / 6
            + (x - 2) ** 4 / 24,
        ),
        ScalarCase(
            "Linear at x=-1",
            "first",
            F=2 * x,
            x0=-1,
            ics_first=0,
            order=3,
            expected=-2 * (x + 1) + (x + 1) ** 2,
        ),
    ]

    # Phase 1: Edge cases (only those with explicit expectations)
    P1_edge = [
        ScalarCase(
            "Order 0",
            "first",
            F=x + sp.Function("y")(x) ** 2,
            x0=0,
            ics_first=5,
            order=0,
            expected=5,
        ),
        ScalarCase("Order 1", "first", F=3 * x, x0=0, ics_first=2, order=1, expected=2),
        ScalarCase(
            "Constant RHS", "first", F=5, x0=0, ics_first=1, order=3, expected=1 + 5 * x
        ),
        ScalarCase(
            "Complex coefficient",
            "first",
            F=I * sp.Function("y")(x),
            x0=0,
            ics_first=1,
            order=3,
            expected=1 + I * x + (I**2) * x**2 / 2 + (I**3) * x**3 / 6,
        ),
    ]

    # Phase 2: Second-order ODEs
    P2_second = [
        ScalarCase(
            "Hyperbolic cosine (y''=y)",
            "nth",
            m=2,
            G=sp.Function("y")(x),
            x0=0,
            ics_nth={0: 1, 1: 0},
            order=6,
            expected=1 + x**2 / 2 + x**4 / 24 + x**6 / 720,
        ),
        ScalarCase(
            "Trigonometric cosine (y''=-y)",
            "nth",
            m=2,
            G=-sp.Function("y")(x),
            x0=0,
            ics_nth={0: 1, 1: 0},
            order=6,
            expected=1 - x**2 / 2 + x**4 / 24 - x**6 / 720,
        ),
        ScalarCase(
            "Trigonometric sine (y''=-y, sine IC)",
            "nth",
            m=2,
            G=-sp.Function("y")(x),
            x0=0,
            ics_nth={0: 0, 1: 1},
            order=5,
            expected=x - x**3 / 6 + x**5 / 120,
        ),
        ScalarCase(
            "SHO frequency 2 (y''=-4y)",
            "nth",
            m=2,
            G=-4 * sp.Function("y")(x),
            x0=0,
            ics_nth={0: 1, 1: 0},
            order=4,
            expected=1 - 2 * x**2 + 2 * x**4 / 3,
        ),
        # Damped oscillator is checked numerically in Maxima - no explicit expected series
    ]

    # Phase 2: Higher-order
    P2_higher = [
        ScalarCase(
            "Third order polynomial (y'''=0)",
            "nth",
            m=3,
            G=0,
            x0=0,
            ics_nth={0: 1, 1: 2, 2: 3},
            order=3,
            expected=1 + 2 * x + (sp.Rational(3, 2)) * x**2,
        ),
        ScalarCase(
            "Fourth order y^(4)=y",
            "nth",
            m=4,
            G=sp.Function("y")(x),
            x0=0,
            ics_nth={0: 1, 1: 0, 2: 0, 3: 0},
            order=8,
            expected=1 + x**4 / 24 + x**8 / 40320,
        ),
        ScalarCase(
            "Third order coupling (y'''=x+y)",
            "nth",
            m=3,
            G=x + sp.Function("y")(x),
            x0=0,
            ics_nth={0: 1, 1: 0, 2: 0},
            order=5,
            expected=1 + x**3 / 6 + x**4 / 24,
        ),
    ]

    # Phase 2: Nonzero expansion for nth-order
    P2_nonzero = [
        ScalarCase(
            "Second order at x=1 (y''=y)",
            "nth",
            m=2,
            G=sp.Function("y")(x),
            x0=1,
            ics_nth={0: 1, 1: 1},
            order=4,
            expected=1
            + (x - 1)
            + (x - 1) ** 2 / 2
            + (x - 1) ** 3 / 6
            + (x - 1) ** 4 / 24,
        ),
        ScalarCase(
            "Constant third derivative at x=2 (y'''=6)",
            "nth",
            m=3,
            G=6,
            x0=2,
            ics_nth={0: 0, 1: 0, 2: 0},
            order=4,
            expected=(x - 2) ** 3,
        ),
    ]

    # Phase 3: Systems with explicit expectations
    f, g = sp.Function("f")(x), sp.Function("g")(x)
    P3_systems = [
        SystemCase(
            "Simple oscillator f/g",
            F=g,
            G=-f,
            x0=0,
            f0=0,
            g0=1,
            order=5,
            expected_f=x - x**3 / 6 + x**5 / 120,
            expected_g=1 - x**2 / 2 + x**4 / 24,
        ),
        SystemCase(
            "Exp system f/g",
            F=f + g,
            G=f + g,
            x0=0,
            f0=1,
            g0=0,
            order=5,
            # NOTE: Maxima test currently expects x^5/15 for both f and g,
            # but the correct coefficient is 2/15.
            expected_f=1
            + x
            + x**2
            + sp.Rational(2, 3) * x**3
            + sp.Rational(1, 3) * x**4
            + sp.Rational(1, 15) * x**5,
            expected_g=x
            + x**2
            + sp.Rational(2, 3) * x**3
            + sp.Rational(1, 3) * x**4
            + sp.Rational(1, 15) * x**5,
        ),
        SystemCase(
            "Mixed poly f/g",
            F=x,
            G=f + g,
            x0=0,
            f0=0,
            g0=1,
            order=4,
            expected_f=x**2 / 2,
            expected_g=1 + x + x**2 / 2 + x**3 / 3 + x**4 / 12,
        ),
        SystemCase(
            "Nonlinear exact (f=0,g=1)",
            F=f * g,
            G=-f,
            x0=0,
            f0=0,
            g0=1,
            order=6,
            expected_f=0,
            expected_g=1,
        ),
    ]

    # Special functions
    P_special = [
        ScalarCase(
            "Integration of cosine",
            "first",
            F=sp.cos(x),
            x0=0,
            ics_first=0,
            order=7,
            expected=x - x**3 / 6 + x**5 / 120 - x**7 / 5040,
        ),
        # "exp" case is numeric check in Maxima, skip explicit expected
        ScalarCase(
            "Integration gives arctan",
            "first",
            F=1 / (1 + x**2),
            x0=0,
            ics_first=0,
            order=9,
            expected=x - x**3 / 3 + x**5 / 5 - x**7 / 7 + x**9 / 9,
        ),
        # The Maxima test labels this as "Integration gives arcsin" but sets F = sqrt(1-x^2).
        # That integrand actually yields: ∫ sqrt(1-x^2) dx = x - x^3/6 - x^5/40 + ...
        ScalarCase(
            "Integration gives arcsin (as written uses sqrt)",
            "first",
            F=sp.sqrt(1 - x**2),
            x0=0,
            ics_first=0,
            order=5,
            expected=x + x**3 / 6 + sp.Rational(3, 40) * x**5,
        ),
    ]

    scalar_cases = (
        P1_basic
        + P1_nonlinear
        + P1_nonzero
        + P1_edge
        + P2_second
        + P2_higher
        + P2_nonzero
        + P_special
    )
    return scalar_cases, P3_systems


# ----------------------------------------------------------------------------
# Runner
# ----------------------------------------------------------------------------


def compare_polys(got: sp.Expr, expect: sp.Expr) -> Tuple[bool, sp.Expr]:
    diff = sp.simplify(sp.expand(got - expect))
    return (diff == 0), diff


def run_scalar_case(c: ScalarCase) -> Tuple[bool, Optional[sp.Expr], Optional[sp.Expr]]:
    # Try exact; otherwise use series recursion
    if c.ode_type == "first":
        # If F is independent of y, direct integration via exact solve should work; but series method is general
        got = series_from_exact(
            sp.Eq(sp.Function("y")(x).diff(x), c.F),
            sp.Function("y")(x),
            {sp.Function("y")(x).subs(x, c.x0): c.ics_first},
            c.order,
            c.x0,
        )
        if got is None:
            got = series_first_order(c.F, c.ics_first, c.order, c.x0)
    else:  # nth
        # rewrite into y'' = G etc in series form directly
        if c.m is None or c.G is None or c.ics_nth is None:
            raise ValueError(f"Malformed nth-order case: {c.label}")
        # exact attempt
        ode_eq = sp.Eq(sp.Function("y")(x).diff(x, c.m), c.G)
        ics_map = {
            sp.Function("y")(x).diff(x, k).subs(x, c.x0): v
            for k, v in c.ics_nth.items()
        }
        got = series_from_exact(ode_eq, sp.Function("y")(x), ics_map, c.order, c.x0)
        if got is None:
            got = series_nth_order(c.m, c.G, c.ics_nth, c.order, c.x0)

    if c.expected is None:
        return True, got, None  # nothing to compare
    ok, diff = compare_polys(got, c.expected)
    return ok, got, diff


def run_system_case(
    sc: SystemCase,
) -> Tuple[bool, Tuple[sp.Expr, sp.Expr], Tuple[Optional[sp.Expr], Optional[sp.Expr]]]:
    # Try exact, else series
    f, g = sp.Function("f")(x), sp.Function("g")(x)
    odes = [sp.Eq(f.diff(x), sc.F), sp.Eq(g.diff(x), sc.G)]
    ics = {f.subs(x, sc.x0): sc.f0, g.subs(x, sc.x0): sc.g0}

    got = series_from_exact(odes, [f, g], ics, sc.order, sc.x0)
    if got is None:
        got = series_system_2x2(sc.F, sc.G, sc.f0, sc.g0, sc.order, sc.x0)  # type: ignore

    # Normalize return
    if isinstance(got, list):
        got_f, got_g = got[0], got[1]
    else:
        got_f, got_g = got  # type: ignore

    ok_f, diff_f = (True, None)
    ok_g, diff_g = (True, None)
    if sc.expected_f is not None:
        ok_f, diff_f = compare_polys(got_f, sc.expected_f)
    if sc.expected_g is not None:
        ok_g, diff_g = compare_polys(got_g, sc.expected_g)

    return (ok_f and ok_g), (got_f, got_g), (diff_f, diff_g)


def main():
    scalar_cases, system_cases = build_cases()

    print("=" * 72)
    print("SYMPY CROSS-CHECK OF MAXIMA TEST EXPECTATIONS")
    print("=" * 72)

    failures: List[str] = []

    print("\n--- Scalar ODE Cases ---")
    for c in scalar_cases:
        ok, got, diff = run_scalar_case(c)
        tag = "PASS" if ok else "FAIL"
        print(f"[{tag}] {c.label}")
        if not ok:
            failures.append(c.label)
            print("  Got     :", poly_to_maxima_str(got))
            print("  Expected:", poly_to_maxima_str(c.expected))
            print("  Diff    :", poly_to_maxima_str(diff))
            # Contextual suggestions
            if "arcsin" in c.label:
                print(
                    "  SUGGESTION: Your test sets y' = sqrt(1 - x^2) but expects arcsin(x) series."
                )
                print(
                    "              Either change the RHS to 1/sqrt(1 - x^2), or update the expected"
                )
                print(
                    "              series to x - x^3/6 - x^5/40 (and rename the test)."
                )
        # Optional: show computed series even when PASS for transparency
        # else:
        #     print("  Series  :", poly_to_maxima_str(got))

    print("\n--- 2×2 System Cases ---")
    for sc in system_cases:
        ok, (got_f, got_g), (diff_f, diff_g) = run_system_case(sc)
        tag = "PASS" if ok else "FAIL"
        print(f"[{tag}] {sc.label}")
        if not ok:
            failures.append(sc.label)
            if diff_f is not None:
                print("  (f) Got     :", poly_to_maxima_str(got_f))
                print("      Expected:", poly_to_maxima_str(sc.expected_f))
                print("      Diff    :", poly_to_maxima_str(diff_f))
            if diff_g is not None:
                print("  (g) Got     :", poly_to_maxima_str(got_g))
                print("      Expected:", poly_to_maxima_str(sc.expected_g))
                print("      Diff    :", poly_to_maxima_str(diff_g))
            if "Exp system" in sc.label:
                print("  SUGGESTION: For f'=f+g, g'=f+g with f(0)=1, g(0)=0,")
                print(
                    "              f=(1+e^{2x})/2, g=(e^{2x}-1)/2 => x^5 coefficient is 2/15 for both,"
                )
                print("              not 1/15. Update the expected series.")

    print("\n" + "-" * 72)
    if failures:
        print("RESULT: SOME TEST ASSUMPTIONS DISAGREE WITH SYMPY")
        print("Offending cases:")
        for name in failures:
            print("  -", name)
        print("\nRecommended Maxima test suite fixes:")
        print(
            '  1) Phase 3 Systems → "Exp system f/g": change x^5 coefficient from 1/15 to 2/15.'
        )
        print('  2) Special Functions → "Integration gives arcsin":')
        print("     EITHER change RHS to 1/sqrt(1-x^2), OR keep RHS sqrt(1-x^2)")
        print(
            "     and change the expected polynomial to x - x^3/6 - x^5/40; also rename the test."
        )
        print(
            "\nPerformance test tip (Maxima): avoid length(string(res)) since 'length' expects a list."
        )
        print("  Use something robust instead, e.g.:")
        print("      length(args(expand(res))) > 10   /* number of summed terms */")
        print("    or hipow(res, x) >= 15             /* highest power check */")
    else:
        print("RESULT: ALL EXPECTATIONS MATCH SYMPY (nice!)")


if __name__ == "__main__":
    main()
