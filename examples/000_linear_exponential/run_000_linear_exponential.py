# Example 000: y' = y, y(0)=1, series about x=0 to order 8

from sympy import symbols, Function, Eq, dsolve, series, exp

x = symbols("x")
y = Function("y")

# Exact solution then series
sol = dsolve(Eq(y(x).diff(x), y(x)), ics={y(0): 1}).rhs
ser = series(sol, x, 0, 9).removeO()  # up to x^8

print("SymPy series (order 8):", ser)
