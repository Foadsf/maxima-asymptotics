# Example 000 — Linear Exponential (y' = y)

We compute the truncated Taylor series for the IVP
\[
y'(x) = y(x),\quad y(0)=1
\]
about \(x=0\) up to order 8. The exact solution is \(y=e^x\), whose series is
\[
1 + x + \frac{x^2}{2} + \frac{x^3}{6} + \frac{x^4}{24} + \frac{x^5}{120}
+ \frac{x^6}{720} + \frac{x^7}{5040} + \frac{x^8}{40320} + \cdots
\]

## How to run

### Maxima
From this folder:
```bat
"C:\maxima-5.48.1\bin\maxima.bat" -b run_000_linear_exponential.mac
````

### Wolfram Mathematica

```bat
"C:\Program Files\Wolfram Research\Wolfram\14.3\wolframscript.exe" -file run_000_linear_exponential.wl
```

### Python (SymPy)

```bat
python run_000_linear_exponential.py
```

All three should print the same polynomial (up to order 8).
