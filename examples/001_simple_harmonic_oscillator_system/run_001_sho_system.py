# Example 001: SHO system f' = g, g' = -f; f(0)=0, g(0)=1; series to order 7

from sympy import symbols, Function, Eq, dsolve, series, sin, cos

x = symbols("x")
f = Function("f")

# Reduce to single ODE: f'' + f = 0, f(0)=0, f'(0)=1
sol_f = dsolve(
    Eq(f(x).diff(x, 2) + f(x), 0), ics={f(0): 0, f(x).diff(x).subs(x, 0): 1}
).rhs  # sin(x)

# g = f'
sol_g = sol_f.diff(x)  # cos(x)

fser = series(sol_f, x, 0, 8).removeO()  # up to x^7
gser = series(sol_g, x, 0, 8).removeO()  # up to x^7

print("SymPy series f (order 7):", fser)
print("SymPy series g (order 7):", gser)
