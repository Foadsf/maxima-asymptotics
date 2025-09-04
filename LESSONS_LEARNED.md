# Lessons Learned — Maxima Asymptotics (Phase 1 MVP)

## 1) Noun vs. verb derivatives
**Issue:** Differentiation and substitution produced nested noun derivatives like
`'diff('diff(y(x),x,1), x, 1)` and even orders like `1+2`, which broke pattern matches.

**Fixes**
- **Canonicalization pass** using `buildq` so orders are numeric:
  - `'diff(y(x),x)` → `'diff(y(x),x,1)`
  - `'diff('diff(y(x),x,n), x, 1)` → `'diff(y(x),x,n+1)` with `n+1` evaluated.
- Multiple passes per step to fully collapse nesting.

## 2) Substitution order (the big one)
**Issue:** Substituting `y'(x)=const` too early created artifacts like `'diff(1,x,1)` and lost higher-order contributions.

**Fix (two-stage substitution per iteration)**
1. **High-order rules first** (for `y''(x)`, `y'''(x)`, …), iterate to a fixed point.
2. Canonicalize again, **then** apply **first-derivative** rules once.
3. Substitute `y(x)=y0`.
4. Zero any derivatives of constants: `'diff(y0,x,n)=0`, `'diff(1,x,n)=0`.
5. Finally set `x=x0`.

Result: correct values like `y'''(0)=8` for the spec example, no stray `'diff(...)`.

## 3) Numeric loop bounds
**Issue:** Bounds like `k+2` sometimes stayed symbolic (e.g., `2+2`), causing “Unable to evaluate predicate …” in loops.

**Fix:** A tiny helper coerces bounds to integers under `simp:true`, ensuring `for` loops always get numeric limits.

## 4) Keep the final series readable
**Issue:** Running `rat/ratsimp` on the entire polynomial recombined terms into a single big fraction.

**Fix:** Only `ratsimp` the **scalar** derivative values; never the whole series. Use `expand` at the end to open `(x-x0)^k`. With `display2d:false`, output is clean and pasteable.

## 5) Avoid array/list pitfalls
**Issue:** `make_array`/`makelist` occasionally received non-atomic integers (e.g., `3+1`) depending on evaluation timing.

**Fix:** Don’t store coefficients in arrays/lists; **accumulate directly** into the solution polynomial.

## 6) Debugging that actually helps
A simple `asym_debug:true$` flag plus `printf` traces (“parse: …”, “k=…: …”) made it obvious where noun diffs or substitutions misfired.

## 7) Validator nuance
**Key point:** For a degree-`n` truncated series, the residual’s derivatives vanish **up to `n-1`**, not `n`. The validator checks 0..`order-1`.

## 8) Negative tests and error surfaces
We used `errcatch(...) = []` to assert that:
- the initial condition matches the expansion point,
- the ODE is solvable for `y'`,
- `order` is a non-negative integer.

These produce expected error messages during the suite but culminate in **ALL PASS ✅**.

## 9) Nonzero expansion points
Always substitute derivative values **before** replacing `y(x)` with `y0` and **before** setting `x=x0`. Then clean up derivatives of constants.

## 10) Portability
Tested on Maxima 5.48.1 (SBCL 2.5.7, Windows). The canonicalization + two-stage approach is resilient across versions because it reduces reliance on internal diff shapes.
