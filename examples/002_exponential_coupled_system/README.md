# Example 002 — Exponential-Coupled System

We consider the linear system
$$
\begin{cases}
f'(x) = f(x) + g(x),\\
g'(x) = f(x) + g(x),
\end{cases}
\qquad f(0)=1,\; g(0)=0.
$$

This decouples via $h=f+g\Rightarrow h'=2h\Rightarrow h=e^{2x}$ and
$v=f-g\Rightarrow v'=0\Rightarrow v=1$. Hence
$$
f(x)=\frac{1+e^{2x}}{2},\qquad g(x)=\frac{e^{2x}-1}{2}.
$$

The series about $x=0$ (to order 7) are
$$
\begin{aligned}
f(x)&=1+x+x^2+\frac{2}{3}x^3+\frac{1}{3}x^4+\frac{2}{15}x^5+\frac{1}{45}x^6+\frac{4}{315}x^7+\cdots,\\
g(x)&=x+x^2+\frac{2}{3}x^3+\frac{1}{3}x^4+\frac{2}{15}x^5+\frac{1}{45}x^6+\frac{4}{315}x^7+\cdots.
\end{aligned}
$$

## How to run

### Maxima
```bat
"C:\maxima-5.48.1\bin\maxima.bat" -b run_002_exponential_coupled_system.mac
````

### Wolfram Mathematica

```bat
"C:\Program Files\Wolfram Research\Wolfram\14.3\wolframscript.exe" -file run_002_exponential_coupled_system.wl
```

### Python (SymPy)

```bat
python run_002_exponential_coupled_system.py
```
