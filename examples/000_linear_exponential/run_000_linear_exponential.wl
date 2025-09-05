(* Example 000: y' = y, y(0)=1, series about x=0 to order 8 *)

ClearAll["Global`*"];
fmt[expr_] := ToString[expr // Expand // Simplify, InputForm];

x = Symbol["x"];
ser = Normal @ AsymptoticDSolveValue[
  {y'[x] == y[x], y[0] == 1},
  y[x], {x, 0, 8}
];

Print["Mathematica series (order 8): ", fmt[ser]];
