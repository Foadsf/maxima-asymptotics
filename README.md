# Maxima Asymptotics Library — Taylor-Series ODE Solver (Phase 1–3)

A compact Maxima library for computing **truncated power-series solutions** of ODEs around a given expansion point, with clean, copy-pasteable polynomials.

* ✅ **Phase 1:** Explicit **first-order** ODEs
* ✅ **Phase 2:** Explicit **nth-order** ODEs with **n** initial conditions
* ✅ **Phase 3:** **Systems of first-order** ODEs (any number of dependent functions)

> **Tested on:** Maxima 5.48.1 (SBCL 2.5.7) on Windows
> **Main file:** `asymptotics.mac`
> **Full test suite:** `comprehensive_test_suite.mac`

---

## Highlights

* **Core (1st-order):** `asymptotic_ode_solve(...)`
* **Nth-order:** `asymptotic_ode_solve_nth(...)`
* **Systems (NEW):** `asymptotic_system_solve(...)`
* **Validation helpers:**

  * `asymptotic_ode_check(...)` (1st-order)
  * `asymptotic_ode_check_nth(...)` (nth-order)
  * `asymptotic_ode_system_check(...)` (systems, NEW)
* **Robust symbolic handling:** canonicalizes noun derivatives and uses **two-stage substitution** to avoid artifacts like `'diff(1,x,1)`.
* **Human-friendly output:** a true sum with rational coefficients (no “one big fraction”). Works great with `display2d:false`.
* **Debugging:** set `asym_debug:true$` to print the internal steps.
* **Cross-checks:** SymPy & Mathematica verifiers in `dev/verification/`.
* **Examples:** runnable, language-triangulated examples in `examples/` (Maxima + Mathematica + SymPy).

---

## Install

1. Put `asymptotics.mac` in your working directory (or on Maxima’s search path).
2. In Maxima:

   ```maxima
   load("asymptotics.mac")$
   display2d:false$     /* optional: linear output for easy copy/paste */
   /* optional debug */
   /* asym_debug:true$ */
   ```

---

## Quickstart

### First-order IVP (about x = 0)

```maxima
asymptotic_ode_solve('diff(y(x),x) = x + y(x)^2, y(x), x=0, 3, [y(0)=1]);
/* -> 1 + x + (3/2)*x^2 + (4/3)*x^3 */
```

### Nth-order example: y″ = −y (cosine), y(0)=1, y′(0)=0

```maxima
asymptotic_ode_solve_nth('diff(y(x),x,2) = -y(x), y(x), x=0, 6,
  [y(0)=1, 'diff(y(x),x,1)=0]);
/* -> 1 - x^2/2 + x^4/24 - x^6/720 */
```

### System example: f′ = g, g′ = −f (sin/cos)

```maxima
asymptotic_system_solve(
  ['diff(f(x),x) = g(x), 'diff(g(x),x) = -f(x)],
  [f(x), g(x)], x=0, 7, [f(0)=0, g(0)=1]);
/* -> [f-series, g-series] = [x - x^3/6 + x^5/120 - x^7/5040, 1 - x^2/2 + x^4/24 - x^6/720] */
```

---

## API Reference

### 1) First-order ODE

```maxima
asymptotic_ode_solve(
  equation,              /* explicit 1st-order: 'diff(y(x),x) = H(x,y)  OR  expr = 0 */
  dep_var(indep_var),    /* e.g. y(x) */
  expansion_point,       /* e.g. x = 0 or x = %pi/4 */
  order,                 /* highest power in the series (non-negative integer) */
  initial_condition      /* list with a single equation, e.g. [y(0) = 1] */
)  /* -> polynomial in (x - x0) */
```

**Validator:**

```maxima
asymptotic_ode_check(
  equation, dep_var(indep_var), expansion_point, order, initial_condition
)  /* -> true/false */

```

Checks that the **residual** and its derivatives vanish at the expansion point **through `order - 1`**.

---

### 2) Nth-order ODE

```maxima
asymptotic_ode_solve_nth(
  equation,              /* y^(n) = H(x, y, y', …, y^(n-1)) OR expr = 0 (canonicalized) */
  dep_var(indep_var),    /* e.g. y(x) */
  expansion_point,       /* e.g. x = 0 */
  order,                 /* series degree (highest power) */
  initial_conditions     /* list of n equations at x0: [y(x0)=a0, 'diff(y(x),x,1)=a1, ..., 'diff(y(x),x,n-1)=a_{n-1}] */
)  /* -> polynomial in (x - x0) */
```

**Validator:**

```maxima
asymptotic_ode_check_nth(
  equation, dep_var(indep_var), expansion_point, order, initial_conditions
)  /* -> true/false */
```

Verifies that `r(x) = y_series^(n) - H(x, y_series, …, y_series^(n-1))` and its derivatives vanish at `x0` through **`max(order - n, 0)`**.

---

### 3) System of first-order ODEs (NEW)

```maxima
asymptotic_system_solve(
  equations_list,        /* e.g. ['diff(f(x),x) = F(x,f,g), 'diff(g(x),x) = G(x,f,g)] */
  dep_vars_list,         /* [f(x), g(x), ...] */
  expansion_point,       /* e.g. x = 0 */
  order,                 /* highest power in the series */
  initial_conditions     /* one IC per dependent function, e.g. [f(0)=..., g(0)=...] */
)  /* -> list of polynomials [series_f, series_g, ...] */
```

**Validator:**

```maxima
asymptotic_ode_system_check(
  equations_list, dep_vars_list, expansion_point, order, initial_conditions
)  /* -> true/false */
```

For each component $i$, it checks that

$$
r_i(x) = \frac{d}{dx}\,series_i(x) \;-\; H_i\bigl(x,\; series\_vector(x)\bigr)
$$

and its derivatives vanish at $x_0$ through **`order - 1`**.
The same canonicalization + derivative-isolation pipeline is used for consistency.

---

## Full Test Suite

The canonical test runner lives in `comprehensive_test_suite.mac` and exercises **Phase 1–3** (including positive cases, numerical cross-checks, and negative/error handling).

### Easiest way (batch file)

A small wrapper is included:

```text
run.mac
```

Run it:

```bat
"C:\maxima-5.48.1\bin\maxima.bat" -b run.mac
```

Or in Maxima:

```maxima
load("asymptotics.mac")$
load("comprehensive_test_suite.mac")$
display2d:false$
asym_debug:false$    /* or true for detailed traces */
run_comprehensive_tests()$
```

You should see a summary ending with **ALL TESTS PASSED** when everything is correct.

---

## Examples

Self-contained, language-triangulated examples live under `examples/`. Each folder has:

* `README.md` — problem statement and how to run
* `run_*.mac` — Maxima script using `../../asymptotics.mac`
* `run_*.wl` — Mathematica script (via `wolframscript`)
* `run_*.py` — SymPy script

Currently available:

* `000_linear_exponential` — $y' = y,\; y(0)=1$
* `001_simple_harmonic_oscillator_system` — $f' = g,\; g' = -f$, $f(0)=0,\; g(0)=1$
* `002_exponential_coupled_system` — $f' = f+g,\; g' = f+g$, $f(0)=1,\; g(0)=0$
* `003_integration_gives_arcsin` — $y' = 1/\sqrt{1-x^2},\; y(0)=0 \Rightarrow y=\arcsin x$

> Tip: you can add a convenience runner (e.g., `examples/run_all_examples.bat`) to execute all examples across Maxima/Mathematica/SymPy.

---

## Cross-Verification (Optional, but recommended)

Independent verifiers live in `dev/verification/`:

* `sympy_crosscheck_all.py` — SymPy cross-check of **all** test assumptions
* `mathematica_crosscheck_all.wl` — Mathematica (AsymptoticDSolveValue/DSolve) cross-check
* `dev/verification/run_all_verifiers.bat` — convenience launcher

Run them (from repo root):

```bat
python dev\verification\sympy_crosscheck_all.py
"C:\Program Files\Wolfram Research\Wolfram\14.3\wolframscript.exe" -file dev\verification\mathematica_crosscheck_all.wl
```

These confirm the expected polynomials in the Maxima tests; any mismatch prints a PASS/FAIL with a diff and suggested fix.

---

## Design Notes

* **Derivative handling (all modes):**

  * Canonicalize *noun* derivatives with evaluated orders across **all** dependent functions before solving for first derivatives (systems) or highest derivative (nth-order).
  * Use **two-stage substitution**:

    1. Apply high-order rules repeatedly (to a small fixed point) to simplify nested derivatives and eliminate artifacts (e.g., `'diff(1,x,1)`).
    2. Canonicalize once more, then apply low-order rules.
  * Only after the above: substitute `x = x0`, kill derivatives of constants/IC values up to the working order, and substitute `y(x)=y0`, `f(x)=f0`, etc.
* **Coefficient recursion:**

  * For 1st-order and systems, compute Taylor coefficients order-by-order by matching the series of the RHS at each step.
  * For nth-order, equate the series of $y^{(n)}$ with the series of $H(x,y,\dots)$, using the factorial factor $\frac{(n+k)!}{k!}$ for coefficient mapping.
* **Stability margin:** internally use a small truncation margin (e.g. `mo := max(order+8, 16)`) during matching/canonicalization to avoid accidental term loss.
* **Formatting:** expand the final polynomial; keep rational terms separate (don’t `ratsimp` the whole sum), but `ratsimp` scalar derivative values while evaluating coefficients.
* **Residual checkers:** always **differentiate before substitution** when evaluating residual derivatives at the expansion point.

---

## Limitations

* **Explicit** ODEs only (no implicit forms yet).
* Expansion about **regular points** (Frobenius/singular points not implemented yet).
* PDEs not supported (planned later).
* For systems, RHS must be explicit in the dependent functions; tested thoroughly for 2D, designed for general size.

---

## Troubleshooting

* Turn on `display2d:false$` for clean linear output.
* Set `asym_debug:true$` to see internal steps and rule applications.
* Negative tests in the suite intentionally print error messages; that’s expected.
* If `'diff(...)` leftovers appear, increase the internal margin (see code) or inspect the debug trace to ensure all derivative forms are canonicalized.

---

## Roadmap

* **Phase 4:** Simple 2D PDEs with multivariate series.
* **Phase 5:** Frobenius / regular singular points; boundary conditions; spectral bases.

---

## License

MIT — see `LICENSE`.

---

## Acknowledgements

Thanks to the cross-verification harnesses (SymPy + Mathematica) for keeping the assumptions honest and catching subtle expectation bugs in the test suite.
