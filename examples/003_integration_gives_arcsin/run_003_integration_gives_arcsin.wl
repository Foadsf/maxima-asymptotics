(* Example 003: y' = 1/Sqrt[1 - x^2], y(0)=0; series to order 7 *)

ClearAll["Global`*"];
fmt[expr_] := ToString[expr // Expand // Simplify, InputForm];
x = Symbol["x"];

ser = Normal @ AsymptoticDSolveValue[
  {y'[x] == 1/Sqrt[1 - x^2], y[0] == 0},
  y[x], {x, 0, 7}
];

Print["Mathematica series (order 7): ", fmt[ser]];
