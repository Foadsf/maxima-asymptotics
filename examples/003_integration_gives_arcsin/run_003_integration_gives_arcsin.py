# Example 003: y' = 1/sqrt(1 - x**2), y(0)=0; series to order 7

from sympy import symbols, Function, Eq, dsolve, series, asin, sqrt

x = symbols("x")
y = Function("y")

# dsolve finds y = asin(x) with IC y(0)=0
sol = dsolve(Eq(y(x).diff(x), 1 / sqrt(1 - x**2)), ics={y(0): 0}).rhs
ser = series(sol, x, 0, 8).removeO()  # up to x^7

print("SymPy series (order 7):", ser)
