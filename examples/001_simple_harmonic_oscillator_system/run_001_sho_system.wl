(* Example 001: SHO system f' = g, g' = -f; f(0)=0, g(0)=1; series to order 7 *)

ClearAll["Global`*"];
fmt[expr_] := ToString[expr // Expand // Simplify, InputForm];

x = Symbol["x"];
{fser, gser} = Normal @ AsymptoticDSolveValue[
  {f'[x] == g[x], g'[x] == -f[x], f[0] == 0, g[0] == 1},
  {f[x], g[x]}, {x, 0, 7}
];

Print["Mathematica series f (order 7): ", fmt[fser]];
Print["Mathematica series g (order 7): ", fmt[gser]];
