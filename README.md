# Maxima Asymptotics Library (Phase 1–2)

A Taylor–series ODE solver library for Maxima. It computes truncated power-series solutions for explicit ODEs around a given point, returning clean, copy-pasteable polynomials.
Phase 2 adds support for **explicit nth-order ODEs** with **n initial conditions**.

> **Tested on** Maxima 5.48.1 (SBCL 2.5.7, Windows).
> **File:** `asymptotics.mac`

## Highlights

- **Core solver:** `asymptotic_ode_solve(...)`
- **Nth-order solver (NEW):** `asymptotic_ode_solve_nth(...)`
- **Robust symbolic handling:** canonicalizes *noun* derivatives and applies two-stage substitutions to avoid artifacts like `'diff(1,x,1)`.
- **Human-friendly output:** returns a sum with rational coefficients (no “one big fraction”). Works great with `display2d:false`.
- **Debugging:** toggle `asym_debug:true$` for step-by-step traces.
- **Validation:** `asymptotic_ode_check(...)` verifies residuals up to `order-1`.
- **Test suite:** `asym_mvp_test_suite()` runs a battery of positive/negative tests.
- **Nth-order validation (NEW):** `asymptotic_ode_check_nth(...)` checks the residual up to `order - n`.
- **Phase 2 tests (NEW):** `asym_phase2_test_suite()` exercises cosh/cos/cubic and residuals.

## Install

1. Place `asymptotics.mac` in your Maxima working directory.
2. In Maxima:
   ```maxima
   load("asymptotics.mac")$
   ```

## API

```maxima
asymptotic_ode_solve(equation, dep_var(indep_var), expansion_point, order, initial_condition)
```

* `equation` — explicit 1st-order ODE (e.g. `'diff(y(x),x) = x + y(x)^2` or `... = 0` form).
* `dep_var(indep_var)` — e.g. `y(x)`.
* `expansion_point` — equation, e.g. `x = 0` or `x = %pi/4`.
* `order` — non-negative integer (highest power in the series).
* `initial_condition` — single-element list, e.g. `[y(0) = 1]`.

**Returns:** truncated series as a symbolic polynomial.

### Validator

```maxima
asymptotic_ode_check(equation, dep_var(indep_var), expansion_point, order, initial_condition)
```

Checks that the residual and its derivatives vanish at the expansion point *through `order - 1`*.

### Test suite

```maxima
asym_mvp_test_suite()
```

Runs standard examples and negative tests. Expected noisy error messages are printed during negative tests; the suite ends with **ALL PASS ✅** when everything is OK.

## API (nth-order)

```maxima
asymptotic_ode_solve_nth(
  equation,          /* explicit nth-order ODE: y^(n) = H(x,y,y',...,y^(n-1)) OR expr=0 */
  dep_var(indep_var),/* e.g. y(x) */
  expansion_point,   /* e.g. x = 0 */
  order,             /* series degree (highest power) */
  initial_conditions /* list of n equations: [y(x0)=a0, 'diff(y(x),x,1)=a1, ..., 'diff(y(x),x,n-1)=a_{n-1}] */
)
```

**Notes**

* The solver **detects** the highest derivative order `n` from the equation (canonicalized), then solves for `y^(n)`.
* Initial conditions are interpreted at the expansion point `x0`.
* Output is an expanded polynomial with rational coefficients (friendly for `display2d:false`).

### Validator (nth-order)

```maxima
asymptotic_ode_check_nth(equation, dep_var(indep_var), expansion_point, order, initial_conditions)
```

Checks that `r(x) = y_series^(n) - H(x, y_series, ..., y_series^(n-1))` and its derivatives vanish at `x0` through **`max(order - n, 0)`**.

## Usage

```maxima
load("asymptotics.mac")$
display2d:false$  /* optional: linear output */

asymptotic_ode_solve('diff(y(x),x) = x + y(x)^2, y(x), x=0, 3, [y(0)=1]);
/* 1 + x + (3/2)*x^2 + (4/3)*x^3 */

asymptotic_ode_solve('diff(y(x),x) = y(x), y(x), x=0, 4, [y(0)=1]);
/* 1 + x + x^2/2 + x^3/6 + x^4/24 */

asymptotic_ode_solve('diff(y(x),x) = (x-2)^3, y(x), x=2, 6, [y(2)=3]);
/* x^4/4 - 2*x^3 + 6*x^2 - 8*x + 7 */

/\* y'' =  y, y(0)=1, y'(0)=0  →  cosh x (through x^6) */
asymptotic\_ode\_solve\_nth('diff(y(x),x,2) = y(x), y(x), x=0, 6,
\[y(0)=1, 'diff(y(x),x,1)=0]);
/* 1 + x^2/2 + x^4/24 + x^6/720 \*/

/\* y'' = -y, y(0)=1, y'(0)=0  →  cos x (through x^6) */
asymptotic\_ode\_solve\_nth('diff(y(x),x,2) = -y(x), y(x), x=0, 6,
\[y(0)=1, 'diff(y(x),x,1)=0]);
/* 1 - x^2/2 + x^4/24 - x^6/720 \*/

/\* y''' = 0, y(0)=1, y'(0)=2, y''(0)=3  →  quadratic */
asymptotic\_ode\_solve\_nth('diff(y(x),x,3) = 0, y(x), x=0, 3,
\[y(0)=1, 'diff(y(x),x,1)=2, 'diff(y(x),x,2)=3]);
/* 1 + 2\*x + (3/2)\*x^2 \*/
```

Validate:

```maxima
asymptotic_ode_check('diff(y(x),x) = y(x) - y(x)^2, y(x), x=0, 5, [y(0)=1/2]);
/* true */
```

Run all tests:

```maxima
asym_mvp_test_suite()$
asym_phase2_test_suite()$
```

## Design Notes

* **Derivative handling:** Normalize nested noun derivatives with evaluated orders, then apply substitutions in two stages: (A) all high-order rules to a fixed point, (B) first-derivative rules once. Only afterward substitute `y(x)=y0`, kill derivatives of constants, and set `x=x0`.
* **Formatting:** Only `expand` the final series; never `ratsimp` the whole polynomial (keeps separate rational terms). We `ratsimp` only scalar derivative values during coefficient evaluation.
* **No arrays/lists for coeffs:** terms are accumulated directly, avoiding `makelist`/`make_array` pitfalls.

## Limitations (Phase 1–2)

* Explicit ODEs only (no implicit forms).
* nth-order ODEs: require **n** initial conditions at the expansion point.
* Regular (non-singular) points only (no Frobenius yet).
* No systems of ODEs, no PDEs (see roadmap).

## Roadmap

* **Phase 3:** Systems of first-order ODEs.
* **Phase 4:** Simple 2D PDEs with multivariate series.
* **Phase 5:** Frobenius/singular points, boundary conditions, spectral bases.

## Troubleshooting

* Set `display2d:false$` if you want easy copy/paste output.
* Negative tests in the suite print expected error messages but still pass overall.
* If you see unexpected `'diff(...)` leftovers, enable `asym_debug:true$` and share the trace.

## License

MIT — see `LICENSE`.
