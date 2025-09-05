# Example 002: f' = f + g, g' = f + g; f(0)=1, g(0)=0; series to order 7

from sympy import symbols, Function, Eq, dsolve, exp, series

x = symbols("x")
f = Function("f")
g = Function("g")

# Solve exactly using change of variables:
# f = (1 + exp(2x))/2, g = (exp(2x) - 1)/2
f_exact = (1 + exp(2 * x)) / 2
g_exact = (exp(2 * x) - 1) / 2

fser = series(f_exact, x, 0, 8).removeO()  # up to x^7
gser = series(g_exact, x, 0, 8).removeO()

print("SymPy series f (order 7):", fser)
print("SymPy series g (order 7):", gser)
