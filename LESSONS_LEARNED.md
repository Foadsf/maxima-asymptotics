# Lessons Learned: Maxima Power Series Library Development

This document captures the key insights, mistakes, and solutions encountered during the development of the Maxima Asymptotics Library. These lessons serve as a reference to avoid repeating errors and to guide future development phases.

## Critical Implementation Mistakes and Solutions

### 1. **Variable Substitution Order Error**

**❌ Mistake**: Attempting to differentiate expressions after substituting numerical values for variables.

**🔍 Root Cause**:
```maxima
current_deriv: subst(indep_var = x0, current_deriv),  # Substituted x = 0 first
current_deriv: diff(current_deriv, indep_var),        # Then tried diff w.r.t. 0 → ERROR
```

**Error Message**: `diff: variable must not be a number; found: 0`

**✅ Solution**: Always differentiate symbolically BEFORE any numerical substitution:
```maxima
current_deriv: diff(current_deriv, indep_var),        # Differentiate first
current_deriv: subst(indep_var = x0, current_deriv),  # Then substitute values
```

**🔑 Principle**: Maxima requires symbolic variables for differentiation. Never substitute numerical values before differentiation operations.

### 2. **Incomplete Function Value Substitution**

**❌ Mistake**: Partial substitution leaving expressions like `y(0)` instead of actual numerical values.

**🔍 Root Cause**:
- Used `subst(dep_var_func = y0, expr)` but this only substituted `y(x)`
- Didn't handle `y(x0)` patterns that appeared after point substitution

**✅ Solution**: Multiple targeted substitutions:
```maxima
eval_result: subst(indep_var = x0, eval_result),          # x → x0 first
eval_result: subst(dep_var_func = y0, eval_result),       # y(x) → y0
eval_result: subst(dep_var(x0) = y0, eval_result),        # y(x0) → y0 explicitly
```

**🔑 Principle**: Maxima's pattern matching in `subst()` is literal. Plan substitution sequences carefully to catch all patterns.

### 3. **Chain Rule Implementation Complexity**

**❌ Mistake**: Over-complicating the chain rule for composite function differentiation.

**🔍 Root Cause**: Attempted to manually implement total derivatives with complex symbolic manipulations that were error-prone.

**✅ Solution**: Use Maxima's built-in differentiation with strategic expression building:
```maxima
# Let Maxima handle the chain rule naturally:
current_expr: diff(current_expr, indep_var) +
              diff(current_expr, dep_var_func) * h_expr
```

**🔑 Principle**: Leverage Maxima's strengths. Don't reinvent symbolic calculus when the system already handles it well.

### 4. **Coefficient Indexing Confusion**

**❌ Mistake**: Mixed 0-based mathematical notation with Maxima's 1-based list indexing.

**🔍 Root Cause**:
- Mathematical: `a₀, a₁, a₂, ...` (0-based)
- Maxima lists: `coeffs[1], coeffs[2], coeffs[3], ...` (1-based)
- Led to off-by-one errors and confusion

**✅ Solution**: Establish clear, consistent mapping and document it:
```maxima
# coeffs[1] = a₀, coeffs[2] = a₁, coeffs[k+1] = aₖ
coeffs: makelist(0, k, 0, expansion_order),  # Create list for a₀ through aₙ
coeffs[1]: y0,                               # a₀ = y(x₀)
coeffs[k+1]: derivative_values[k+1] / k!,    # aₖ = y⁽ᵏ⁾(x₀) / k!
```

**🔑 Principle**: Document indexing conventions clearly. Use consistent mapping throughout the codebase.

## Algorithm Design Insights

### 5. **Iterative Derivative Building Strategy**

**✅ Key Insight**: Build derivatives iteratively by maintaining a "current expression" that represents the most recent derivative:

```maxima
current_expr: h_expr,                    # Start with y'(x) = h(x,y)
for k: 2 thru order do (
    current_expr: diff(...),             # Build y''(x), y'''(x), etc.
    # Evaluate at expansion point
    # Store coefficient
)
```

**🔑 Principle**: Sequential building is more reliable than trying to compute high-order derivatives directly.

### 6. **Separation of Symbolic and Numerical Operations**

**✅ Key Insight**: Maintain clear separation between symbolic manipulation and numerical evaluation:

1. **Symbolic Phase**: Build all derivative expressions
2. **Evaluation Phase**: Substitute values and compute coefficients
3. **Construction Phase**: Build the final polynomial

**🔑 Principle**: Don't mix symbolic and numerical operations within the same computational step.

## Input Validation Lessons

### 7. **Comprehensive Error Checking**

**✅ Best Practice**: Validate every input parameter with specific, helpful error messages:

```maxima
if not (integerp(expansion_order) and expansion_order >= 0) then
    error("expansion_order must be a non-negative integer"),
```

**🔑 Principle**: Fail fast with clear error messages. Good error handling is as important as the algorithm itself.

### 8. **Consistent Parameter Parsing**

**✅ Best Practice**: Extract and validate components systematically:
- Use `op()` and `args()` consistently for expression decomposition
- Validate structural requirements before accessing components
- Check consistency between related parameters (expansion point vs. initial condition)

## Maxima-Specific Technical Insights

### 9. **Expression Structure Understanding**

**✅ Key Learning**: Maxima expressions have predictable structure:
- `op(expr)` gives the main operator
- `args(expr)` gives the operands as a list
- Equations (`=`) are operators with left/right-hand sides accessible via `lhs()` and `rhs()`

### 10. **Solve Function Behavior**

**✅ Key Learning**: `solve()` returns a list of solutions, even for single solutions:
```maxima
h_expr: solve(equation, 'diff(dep_var_func, indep_var)),  # Returns [solution]
h_expr: rhs(h_expr[1]),                                   # Extract first solution's RHS
```

## Development Process Insights

### 11. **Incremental Testing Strategy**

**✅ Best Practice**: Test each component in isolation:
1. Input parsing and validation
2. Derivative extraction
3. Coefficient calculation for simple cases
4. Full algorithm with known examples

**🔑 Principle**: Build confidence through incremental verification rather than testing the complete system at once.

### 12. **Manual Verification as Ground Truth**

**✅ Best Practice**: Include manual calculations for test cases:
```maxima
/* Manual verification for y' = x + y², y(0) = 1:
 * a[0] = 1
 * a[1] = h(0,1) = 1
 * a[2] = (1 + 2*1*1)/2! = 3/2
 * a[3] = (2*1² + 2*1*3)/3! = 4/3
 * Expected: 1 + x + (3/2)*x² + (4/3)*x³
 */
```

**🔑 Principle**: Manual verification provides confidence and catches algorithmic errors that might not be obvious from code inspection alone.

## Phase 1 Implementation Lessons

### 13. **Noun vs. verb derivatives**

**Issue:** Differentiation and substitution produced nested noun derivatives like
`'diff('diff(y(x),x,1), x, 1)` and even orders like `1+2`, which broke pattern matches.

**Fixes**
- **Canonicalization pass** using `buildq` so orders are numeric:
  - `'diff(y(x),x)` → `'diff(y(x),x,1)`
  - `'diff('diff(y(x),x,n), x, 1)` → `'diff(y(x),x,n+1)` with `n+1` evaluated.
- Multiple passes per step to fully collapse nesting.

### 14. **Substitution order (the big one)**

**Issue:** Substituting `y'(x)=const` too early created artifacts like `'diff(1,x,1)` and lost higher-order contributions.

**Fix (two-stage substitution per iteration)**
1. **High-order rules first** (for `y''(x)`, `y'''(x)`, …), iterate to a fixed point.
2. Canonicalize again, **then** apply **first-derivative** rules once.
3. Substitute `y(x)=y0`.
4. Zero any derivatives of constants: `'diff(y0,x,n)=0`, `'diff(1,x,n)=0`.
5. Finally set `x=x0`.

Result: correct values like `y'''(0)=8` for the spec example, no stray `'diff(...)`.

### 15. **Numeric loop bounds**

**Issue:** Bounds like `k+2` sometimes stayed symbolic (e.g., `2+2`), causing "Unable to evaluate predicate …" in loops.

**Fix:** A tiny helper coerces bounds to integers under `simp:true`, ensuring `for` loops always get numeric limits.

### 16. **Keep the final series readable**

**Issue:** Running `rat/ratsimp` on the entire polynomial recombined terms into a single big fraction.

**Fix:** Only `ratsimp` the **scalar** derivative values; never the whole series. Use `expand` at the end to open `(x-x0)^k`. With `display2d:false`, output is clean and pasteable.

### 17. **Avoid array/list pitfalls**

**Issue:** `make_array`/`makelist` occasionally received non-atomic integers (e.g., `3+1`) depending on evaluation timing.

**Fix:** Don't store coefficients in arrays/lists; **accumulate directly** into the solution polynomial.

### 18. **Debugging that actually helps**

A simple `asym_debug:true$` flag plus `printf` traces ("parse: …", "k=…: …") made it obvious where noun diffs or substitutions misfired.

### 19. **Validator nuance**

**Key point:** For a degree-`n` truncated series, the residual's derivatives vanish **up to `n-1`**, not `n`. The validator checks 0..`order-1`.

### 20. **Negative tests and error surfaces**

We used `errcatch(...) = []` to assert that:
- the initial condition matches the expansion point,
- the ODE is solvable for `y'`,
- `order` is a non-negative integer.

These produce expected error messages during the suite but culminate in **ALL PASS ✅**.

### 21. **Nonzero expansion points**

Always substitute derivative values **before** replacing `y(x)` with `y0` and **before** setting `x=x0`. Then clean up derivatives of constants.

### 22. **Portability**

Tested on Maxima 5.48.1 (SBCL 2.5.7, Windows). The canonicalization + two-stage approach is resilient across versions because it reduces reliance on internal diff shapes.

## Phase 2 (nth-order) Additions

### 23. **Robust highest-order detection**

**Problem:** Structural scans sometimes missed the true highest derivative, causing the solver to target `y'` instead of `y^(n)`.

**Fix:** A 3-stage detector on a *canonicalized* equation:

1. Structural walk (`asym_max_deriv_order`) on LHS/RHS,
2. `freeof` scan for `'diff(y(x),x,k)` from high to low,
3. `solve` probe for `'diff(y(x),x,k)`.
   Whichever hits first picks `n`. We log `eqn canonical LHS/RHS`, the chosen stage, and the final `n`.

### 24. **Canonicalize once, use everywhere**

We canonicalize the equation to normalize noun derivatives and numeric orders, then reuse that canonical form for both detection and the `solve` that extracts `H`. This avoids drift between detection and solve targets.

### 25. **IC parsing: assoc returns the value**

In Maxima, `assoc(k, alist)` returns the value (or `false`), not a key—value pair. Use `ak : assoc(k, ic_vals)` directly. Also keep substitution rules as a **list** (`listify(setify(...))`) so `subst` receives the right type.

### 26. **Residual order (nth-order)**

For a degree-`N` series of an `n`-th order ODE, the residual derivatives vanish at the expansion point through **`max(N - n, 0)`** (first-order is the special case `N - 1`). The new `asymptotic_ode_check_nth(...)` enforces this.

### 27. **Debugging that scales**

We extended `asym_debug` traces with:

* canonical LHS/RHS of the equation,
* detector stage outcomes,
* the exact derivative variable we solve for,
* the evaluated `y^(k)(x0)` per step.
  This made failures immediately localizable.

## Future Development Guidelines

### 28. **Extensibility Considerations**

Based on the current implementation, future phases should:

- **Maintain the coefficient-based approach** - it scales well to higher orders
- **Preserve the input validation pattern** - it's comprehensive and user-friendly
- **Extend the derivative building strategy** - it can handle multiple variables and higher orders
- **Keep symbolic/numerical separation** - it will be crucial for PDEs and systems

### 29. **Code Organization Principles**

- **Single responsibility**: Each function should have one clear purpose
- **Clear variable naming**: `expansion_order` vs. `order`, `dep_var_func` vs. `y`
- **Comprehensive commenting**: Explain both the mathematics and the implementation
- **Error handling first**: Validate inputs before processing


## Phase 3 (Systems) Lessons

### 30. Canonicalize across *all* dependent functions

When working with systems, always canonicalize noun derivatives for **every** dependent symbol before attempting to isolate first derivatives. This prevents mismatched forms like `'diff(f(x),x,1)` coexisting with `'diff('diff(f(x),x,1), x, 1)`. The two-stage substitution strategy from Phase 1 (high-order to fixed point, then low-order once) **must** be applied component-wise to be reliable. (Builds on §14 and §21.)&#x20;

### 31. Isolate *all* first derivatives once, then recurse

For systems $(\mathbf{y}'=\mathbf{H}(x,\mathbf{y}))$, first use `solve` to extract the full RHS vector $\mathbf{H}$ in one shot on the **canonicalized** equations, then do series coefficient recursion jointly. Re-solving per component per order is slower and can drift structurally. (Related to §10 “Solve function behavior”.)&#x20;

### 32. Two-stage rules that actually converge (systems)

The **Stage-A/Stage-B** rule application works for systems too:

1. Stage-A (iterate): collapse nested derivatives and apply **high-order** rules to a fixed point with a small safety margin $mo=\max(N+8,16)$.
2. Stage-B (single pass): canonicalize again, then apply **low-order** / first-derivative rules once.
   Only after that, substitute values for $x=x_0$, function values $f_i(x_0)$, and kill derivatives of constants. (Extends §14 and §21.)&#x20;

### 33. Residual order for systems

For a degree-$N$ system series, each component residual and its derivatives must vanish at $x_0$ through **$N-1$** (just like the scalar first-order case in §19). Our `asymptotic_ode_system_check(...)` implements exactly that.&#x20;

### 34. Input validation that saves you later

The error “equations\_list and dep\_vars\_list must have the same length” is **expected** in negative tests; keep strict checks: equal list lengths, one IC per dependent function, single shared independent variable, non-negative `order`, ICs evaluated at the **same** expansion point, etc. (Expands §7–§8, §20.)&#x20;

---

## Cross-Verification & Assumption Drift

### 35. Independent verifiers prevent “assumption rot”

Using SymPy *and* Mathematica surfaced two important expectation bugs in our test suite:

* **System $f'=f+g,\;g'=f+g$**: the $x^5$ coefficient is **$2/15$** (not $1/15$).
* **“Integration gives arcsin”**: the test used $y'=\sqrt{1-x^2}$ but expected the **arcsin** series (whose derivative is $1/\sqrt{1-x^2}$).
  Cross-checks caught both and kept our Maxima tests honest. (Supports the “Manual Verification” idea in §12 and expands the general testing story.)&#x20;

### 36. Mathematica quirks in CLI mode

`AsymptoticDSolveValue` sometimes fails for systems under `wolframscript` even when `DSolveValue` works. A robust fallback is a **manual coefficient recursion** for 2×2 systems (what we implemented in our WL harness). Prefer that path first for deterministic CLI runs. (New lesson.)

### 37. Performance checks: don’t `length(string(...))`

Maxima’s `length` expects a **list**, not a string. Use list/structure metrics instead:

* `length(args(expand(res)))  /* term count */`
* `hipow(res, x)              /* highest power */`
  This avoids type errors in performance tests. (New lesson.)

---

## Examples & Repository Hygiene

### 38. Triangulated examples increase trust

Each example folder now includes a **Maxima** script, a **SymPy** script, and a **Mathematica** script that compute the *same* series. This “triangle” makes it much easier to diagnose discrepancies—if one corner disagrees, you know where to look. (New lesson.)

### 39. Keep verifiers out of the root

Move cross-check tools into `dev/verification/` (with a `legacy/` corner). The root should contain the library, the suite, and minimal run files. A tiny `run_all_verifiers.bat` makes it trivial to re-validate everything before commits. (New lesson.)

---

## Concrete Snippets to Add (for quick copy/paste)

**Systems residual principle** (documentation reminder):

```maxima
/* For a degree-N system series, verify residuals r_i and r_i^(j)(x0) = 0
   for all components i and all j = 0..(N-1). */
```

**Performance metric** (replace string-length checks):

```maxima
/* Prefer structure-based checks, not strings */
length(args(expand(res))) > 10  /* number of terms */
hipow(res, x) >= 15             /* highest power */
```

**Mathematica fallback hint** (comment in WL harness):

```wl
(* For systems, prefer manual series recursion under wolframscript to avoid
   sporadic AsymptoticDSolveValue/DSolveValue failures in CLI. *)
```

## Summary

The key to successful Maxima library development is understanding the system's symbolic nature and working with it, not against it. The most critical lesson is **differentiate first, substitute later** - this single principle could have prevented the majority of debugging time spent on this project.

The two-stage substitution pattern (high-order rules to fixed point, then low-order rules once) emerged as the most reliable approach for handling complex derivative expressions without artifacts.

These lessons should guide all future development phases and help maintain the quality and reliability of the asymptotics library as it grows in scope and complexity.
