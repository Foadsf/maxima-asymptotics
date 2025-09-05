# Example 001 — Simple Harmonic Oscillator (2×2 system)

We study the first-order system
$$
\begin{cases}
f'(x) = g(x),\\
g'(x) = -f(x),
\end{cases}
\quad f(0)=0,\ g(0)=1,
$$
whose exact solution is $f(x)=\sin x,\ g(x)=\cos x$.

About $x=0$, to order 7:
$$
f(x) = x - \frac{x^3}{6} + \frac{x^5}{120} - \frac{x^7}{5040} + \cdots
$$
$$
g(x) = 1 - \frac{x^2}{2} + \frac{x^4}{24} - \frac{x^6}{720} + \cdots
$$

## How to run

### Maxima
```bat
"C:\maxima-5.48.1\bin\maxima.bat" -b run_001_sho_system.mac
````

### Wolfram Mathematica

```bat
"C:\Program Files\Wolfram Research\Wolfram\14.3\wolframscript.exe" -file run_001_sho_system.wl
```

### Python (SymPy)

```bat
python run_001_sho_system.py
```

All three should produce matching series up to the requested order.
