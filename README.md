# Maxima Asymptotics Library (Phase 1 MVP)

A Taylor–series ODE solver library for Maxima. It computes truncated power-series solutions for **single, first-order, explicit ODEs** with an initial condition, returning a clean, copy-pasteable polynomial.

> **Tested on** Maxima 5.48.1 (SBCL 2.5.7, Windows).
> **File:** `asymptotics.mac`

## Highlights

- **Core solver:** `asymptotic_ode_solve(...)`
- **Robust symbolic handling:** canonicalizes *noun* derivatives and applies two-stage substitutions to avoid artifacts like `'diff(1,x,1)`.
- **Human-friendly output:** returns a sum with rational coefficients (no “one big fraction”). Works great with `display2d:false`.
- **Debugging:** toggle `asym_debug:true$` for step-by-step traces.
- **Validation:** `asymptotic_ode_check(...)` verifies residuals up to `order-1`.
- **Test suite:** `asym_mvp_test_suite()` runs a battery of positive/negative tests.

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
```

Validate:

```maxima
asymptotic_ode_check('diff(y(x),x) = y(x) - y(x)^2, y(x), x=0, 5, [y(0)=1/2]);
/* true */
```

Run all tests:

```maxima
asym_mvp_test_suite()$
```

## Design Notes

* **Derivative handling:** Normalize nested noun derivatives with evaluated orders, then apply substitutions in two stages: (A) all high-order rules to a fixed point, (B) first-derivative rules once. Only afterward substitute `y(x)=y0`, kill derivatives of constants, and set `x=x0`.
* **Formatting:** Only `expand` the final series; never `ratsimp` the whole polynomial (keeps separate rational terms). We `ratsimp` only scalar derivative values during coefficient evaluation.
* **No arrays/lists for coeffs:** terms are accumulated directly, avoiding `makelist`/`make_array` pitfalls.

## Limitations (MVP)

* Single ODE, first order, explicit, one initial condition.
* Regular (non-singular) points only.
* No systems, higher-order ODEs, or PDEs (see roadmap).

## Roadmap

* **Phase 2:** nth-order explicit ODEs (with n initial conditions).
* **Phase 3:** systems of first-order ODEs.
* **Phase 4:** simple 2D PDEs (multivariate series).
* **Phase 5:** Frobenius/singular points, boundary conditions, spectral bases.

## Troubleshooting

* Set `display2d:false$` if you want easy copy/paste output.
* Negative tests in the suite print expected error messages but still pass overall.
* If you see unexpected `'diff(...)` leftovers, enable `asym_debug:true$` and share the trace.

## License

MIT — see `LICENSE`.
