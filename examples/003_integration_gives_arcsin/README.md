# Example 003 — Integration gives arcsin

We compute the series for the IVP

$$
y'(x) = \frac{1}{\sqrt{1 - x^2}},\quad y(0)=0
$$

so that $y(x) = \arcsin(x)$.

About $x=0$, to order 7:

$$
\arcsin(x) = x + \frac{x^3}{6} + \frac{3x^5}{40} + \frac{5x^7}{112} + \cdots
$$

## How to run

### Maxima
```bat
"C:\maxima-5.48.1\bin\maxima.bat" -b run_003_integration_gives_arcsin.mac
````

### Wolfram Mathematica

```bat
"C:\Program Files\Wolfram Research\Wolfram\14.3\wolframscript.exe" -file run_003_integration_gives_arcsin.wl
```

### Python (SymPy)

```bat
python run_003_integration_gives_arcsin.py
```
