(* ::Package:: *)

(*  sympy_crosscheck_all_mathematica.wl
    -----------------------------------
    Cross-check the *assumptions* (expected Taylor polynomials) in your Maxima
    test suite using Wolfram Mathematica's AsymptoticDSolveValue.

    Run from Windows CMD:
      "C:\Program Files\Wolfram Research\Wolfram\14.3\wolframscript.exe" -file sympy_crosscheck_all_mathematica.wl

    What this script does:
      • Defines the same ODE test cases (scalar & 2×2 systems) with initial conditions
        and expansion orders as in your Maxima/SymPy runs.
      • Uses AsymptoticDSolveValue to compute the series about the specified expansion point.
      • Compares to the hardcoded "expected" polynomials.
      • Prints PASS/FAIL with Got/Expected/Diff, mirroring your Python run.
      • Highlights the two problematic assumptions:
           - Phase 3 "Exp system f/g" (x^5 coefficient should be 2/15)
           - "Integration gives arcsin" (test uses sqrt(1-x^2) but expects arcsin series)

    Notes:
      • All arithmetic uses rationals for exactness.
      • Output uses InputForm for clean CLI printing.
*)

ClearAll["Global`*"];
$Assumptions = True;

fmt[expr_] := ToString[expr // Expand // Simplify, InputForm];

printHeader[title_] := (
  Print[StringRepeat["=", 72]];
  Print[title];
  Print[StringRepeat["=", 72]];
);

printSection[title_] := Print["\n--- ", title, " ---"];

printCase[label_, ok_, got_, expected_, opts:OptionsPattern[]] := Module[{},
  Print["[", If[TrueQ[ok], "PASS", "FAIL"], "] ", label];
  If[TrueQ[ok], Return[]];
  Print["  Got     : ", fmt[got]];
  Print["  Expected: ", fmt[expected]];
  Print["  Diff    : ", fmt[Expand[got - expected]]];
  If[TrueQ[OptionValue["ArcsinHint"]],
    Print["  SUGGESTION: Your test sets y'[x] = Sqrt[1 - x^2] but expects arcsin(x) series."];
    Print["              Either change RHS to 1/Sqrt[1 - x^2], or update the expected"];
    Print["              series to x - x^3/6 - x^5/40 (and rename the test)."];
  ];
  If[TrueQ[OptionValue["ExpSystemHint"]],
    Print["  SUGGESTION: For f' = f + g, g' = f + g with f(0) = 1, g(0) = 0,"];
    Print["              f = (1 + E^(2 x))/2, g = (E^(2 x) - 1)/2 so x^5 coefficient is 2/15 for both,"];
    Print["              not 1/15. Update the expected series."];
  ];
];
Options[printCase] = {"ArcsinHint"->False, "ExpSystemHint"->False};

compare[got_, expected_] := TrueQ[Simplify[Expand[got - expected] === 0]];

(* Robust series solver: AsymptoticDSolveValue, fallbacks if needed *)
scalarSeriesSolve[eqns_List, ics_List, y_, x_, x0_, ord_Integer] := Module[{ser},
  ser = Quiet@Check[
    Normal@AsymptoticDSolveValue[{eqns, ics}, y[x], {x, x0, ord}],
    $Failed
  ];
  If[ser === $Failed,
    ser = Quiet@Check[
      Normal@Series[
        DSolveValue[{eqns, ics}, y[x], x]
        , {x, x0, ord}],
      $Failed
    ];
  ];
  ser
];

systemSeriesSolve2[eqns_List, ics_List, vars_List, x_, x0_, ord_Integer] := Module[{ser},
  ser = Quiet@Check[
    Normal@AsymptoticDSolveValue[{eqns, ics}, vars /. v_Symbol :> v[x], {x, x0, ord}],
    $Failed
  ];
  If[ser === $Failed,
    ser = Quiet@Check[
      Normal /@ (Series[#, {x, x0, ord}] & /@ (DSolveValue[{eqns, ics}, vars /. v_Symbol :> v[x], x])),
      $Failed
    ];
  ];
  ser
];

xSym = x; (* symbol used everywhere *)

(* ----------------------------- *)
(* Build test cases (Scalar ODE) *)
(* ----------------------------- *)
scalarCases = {
  <|
    "Label"->"Linear homogeneous y'=y",
    "Eqns" -> { y'[x] == y[x] },
    "ICs"  -> { y[0]==1 },
    "Var"  -> y, "x0"->0, "Order"->6,
    "Expected" -> 1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 + x^6/720
  |>,
  <|
    "Label"->"Linear y'=2y",
    "Eqns" -> { y'[x] == 2 y[x] },
    "ICs"  -> { y[0]==3 },
    "Var"  -> y, "x0"->0, "Order"->4,
    "Expected" -> 3 + 6 x + 6 x^2 + 4 x^3 + 2 x^4
  |>,
  <|
    "Label"->"Polynomial source y'=x^2",
    "Eqns" -> { y'[x] == x^2 },
    "ICs"  -> { y[0]==0 },
    "Var"  -> y, "x0"->0, "Order"->5,
    "Expected" -> x^3/3
  |>,
  <|
    "Label"->"Mixed polynomial y'=x+x^2",
    "Eqns" -> { y'[x] == x + x^2 },
    "ICs"  -> { y[0]==1 },
    "Var"  -> y, "x0"->0, "Order"->4,
    "Expected" -> 1 + x^2/2 + x^3/3
  |>,
  <|
    "Label"->"Separable y'=y^2",
    "Eqns" -> { y'[x] == y[x]^2 },
    "ICs"  -> { y[0]==1 },
    "Var"  -> y, "x0"->0, "Order"->5,
    "Expected" -> 1 + x + x^2 + x^3 + x^4 + x^5
  |>,
  (* Nonlinear *)
  <|
    "Label"->"Bernoulli y'=xy",
    "Eqns" -> { y'[x] == x y[x] },
    "ICs"  -> { y[0]==1 },
    "Var"  -> y, "x0"->0, "Order"->6,
    "Expected" -> 1 + x^2/2 + x^4/8 + x^6/48
  |>,
  <|
    "Label"->"Linear inhomo y'=y+sin(x)",
    "Eqns" -> { y'[x] == y[x] + Sin[x] },
    "ICs"  -> { y[0]==0 },
    "Var"  -> y, "x0"->0, "Order"->5,
    "Expected" -> x^2/2 + x^3/6
  |>,
  <|
    "Label"->"Riccati y'=1+y^2",
    "Eqns" -> { y'[x] == 1 + y[x]^2 },
    "ICs"  -> { y[0]==0 },
    "Var"  -> y, "x0"->0, "Order"->5,
    "Expected" -> x + x^3/3 + 2 x^5/15
  |>,
  (* Non-zero expansion points *)
  With[{t = (x - 1)},
    <|
      "Label"->"Nonzero point y'=x^2 at x=1",
      "Eqns" -> { y'[x] == x^2 },
      "ICs"  -> { y[1]==2 },
      "Var"  -> y, "x0"->1, "Order"->4,
      "Expected" -> 2 + t + t^2 + t^3/3
    |>
  ],
  <|
    "Label"->"Exponential at x=2",
    "Eqns" -> { y'[x] == y[x] },
    "ICs"  -> { y[2]==1 },
    "Var"  -> y, "x0"->2, "Order"->4,
    "Expected" -> 1 + (x-2) + (x-2)^2/2 + (x-2)^3/6 + (x-2)^4/24
  |>,
  <|
    "Label"->"Linear at x=-1",
    "Eqns" -> { y'[x] == 2 x },
    "ICs"  -> { y[-1]==0 },
    "Var"  -> y, "x0"->-1, "Order"->3,
    "Expected" -> -2 (x+1) + (x+1)^2
  |>,
  (* Edge cases *)
  <|
    "Label"->"Order 0",
    "Eqns" -> { y'[x] == x + y[x]^2 },
    "ICs"  -> { y[0]==5 },
    "Var"  -> y, "x0"->0, "Order"->0,
    "Expected" -> 5
  |>,
  <|
    "Label"->"Order 1",
    "Eqns" -> { y'[x] == 3 x },
    "ICs"  -> { y[0]==2 },
    "Var"  -> y, "x0"->0, "Order"->1,
    "Expected" -> 2
  |>,
  <|
    "Label"->"Constant RHS",
    "Eqns" -> { y'[x] == 5 },
    "ICs"  -> { y[0]==1 },
    "Var"  -> y, "x0"->0, "Order"->3,
    "Expected" -> 1 + 5 x
  |>,
  <|
    "Label"->"Complex coefficient",
    "Eqns" -> { y'[x] == I y[x] },
    "ICs"  -> { y[0]==1 },
    "Var"  -> y, "x0"->0, "Order"->3,
    "Expected" -> 1 + I x + (I^2) x^2/2 + (I^3) x^3/6
  |>,
  (* Second order *)
  <|
    "Label"->"Hyperbolic cosine (y''=y)",
    "Eqns" -> { y''[x] == y[x] },
    "ICs"  -> { y[0]==1, y'[0]==0 },
    "Var"  -> y, "x0"->0, "Order"->6,
    "Expected" -> 1 + x^2/2 + x^4/24 + x^6/720
  |>,
  <|
    "Label"->"Trigonometric cosine (y''=-y)",
    "Eqns" -> { y''[x] == -y[x] },
    "ICs"  -> { y[0]==1, y'[0]==0 },
    "Var"  -> y, "x0"->0, "Order"->6,
    "Expected" -> 1 - x^2/2 + x^4/24 - x^6/720
  |>,
  <|
    "Label"->"Trigonometric sine (y''=-y, sine IC)",
    "Eqns" -> { y''[x] == -y[x] },
    "ICs"  -> { y[0]==0, y'[0]==1 },
    "Var"  -> y, "x0"->0, "Order"->5,
    "Expected" -> x - x^3/6 + x^5/120
  |>,
  <|
    "Label"->"SHO frequency 2 (y''=-4y)",
    "Eqns" -> { y''[x] == -4 y[x] },
    "ICs"  -> { y[0]==1, y'[0]==0 },
    "Var"  -> y, "x0"->0, "Order"->4,
    "Expected" -> 1 - 2 x^2 + (2/3) x^4
  |>,
  (* Higher order *)
  <|
    "Label"->"Third order polynomial (y'''=0)",
    "Eqns" -> { y'''[x] == 0 },
    "ICs"  -> { y[0]==1, y'[0]==2, y''[0]==3 },
    "Var"  -> y, "x0"->0, "Order"->3,
    "Expected" -> 1 + 2 x + (3/2) x^2
  |>,
  <|
    "Label"->"Fourth order y^(4)=y",
    "Eqns" -> { y''''[x] == y[x] },
    "ICs"  -> { y[0]==1, y'[0]==0, y''[0]==0, y'''[0]==0 },
    "Var"  -> y, "x0"->0, "Order"->8,
    "Expected" -> 1 + x^4/24 + x^8/40320
  |>,
  <|
    "Label"->"Third order coupling (y'''=x+y)",
    "Eqns" -> { y'''[x] == x + y[x] },
    "ICs"  -> { y[0]==1, y'[0]==0, y''[0]==0 },
    "Var"  -> y, "x0"->0, "Order"->5,
    "Expected" -> 1 + x^3/6 + x^4/24
  |>,
  (* Nonzero expansion nth-order *)
  <|
    "Label"->"Second order at x=1 (y''=y)",
    "Eqns" -> { y''[x] == y[x] },
    "ICs"  -> { y[1]==1, y'[1]==1 },
    "Var"  -> y, "x0"->1, "Order"->4,
    "Expected" -> 1 + (x-1) + (x-1)^2/2 + (x-1)^3/6 + (x-1)^4/24
  |>,
  <|
    "Label"->"Constant third derivative at x=2 (y'''=6)",
    "Eqns" -> { y'''[x] == 6 },
    "ICs"  -> { y[2]==0, y'[2]==0, y''[2]==0 },
    "Var"  -> y, "x0"->2, "Order"->4,
    "Expected" -> (x-2)^3
  |>,
  (* Special functions *)
  <|
    "Label"->"Integration of cosine",
    "Eqns" -> { y'[x] == Cos[x] },
    "ICs"  -> { y[0]==0 },
    "Var"  -> y, "x0"->0, "Order"->7,
    "Expected" -> x - x^3/6 + x^5/120 - x^7/5040
  |>,
  <|
    "Label"->"Integration gives arctan",
    "Eqns" -> { y'[x] == 1/(1 + x^2) },
    "ICs"  -> { y[0]==0 },
    "Var"  -> y, "x0"->0, "Order"->9,
    "Expected" -> x - x^3/3 + x^5/5 - x^7/7 + x^9/9
  |>,
  (* Intentional mismatch to validate the suite: label says arcsin but RHS is sqrt *)
  <|
    "Label"->"Integration gives arcsin (as written uses sqrt)",
    "Eqns" -> { y'[x] == Sqrt[1 - x^2] },
    "ICs"  -> { y[0]==0 },
    "Var"  -> y, "x0"->0, "Order"->5,
    "Expected" -> x + x^3/6 + 3 x^5/40, (* arcsin series — intentionally wrong for this RHS *)
    "ArcsinHint" -> True
  |>
};

(* ----------------------------- *)
(* Build test cases (2×2 systems) *)
(* ----------------------------- *)
systemCases = {
  <|
    "Label"->"Simple oscillator f/g",
    "Eqns" -> { f'[x] == g[x], g'[x] == -f[x] },
    "ICs"  -> { f[0]==0, g[0]==1 },
    "Vars" -> {f, g}, "x0"->0, "Order"->5,
    "Expected" -> { x - x^3/6 + x^5/120, 1 - x^2/2 + x^4/24 }
  |>,
  <|
    "Label"->"Exp system f/g",
    "Eqns" -> { f'[x] == f[x] + g[x], g'[x] == f[x] + g[x] },
    "ICs"  -> { f[0]==1, g[0]==0 },
    "Vars" -> {f, g}, "x0"->0, "Order"->5,
    "Expected" -> {
      1 + x + x^2 + (2 x^3)/3 + (x^4)/3 + (x^5)/15,  (* <-- hardcoded *wrong* expectation *)
      x + x^2 + (2 x^3)/3 + (x^4)/3 + (x^5)/15       (* correct is 2/15 *)
    },
    "ExpSystemHint" -> True
  |>,
  <|
    "Label"->"Mixed poly f/g",
    "Eqns" -> { f'[x] == x, g'[x] == f[x] + g[x] },
    "ICs"  -> { f[0]==0, g[0]==1 },
    "Vars" -> {f, g}, "x0"->0, "Order"->4,
    "Expected" -> { x^2/2, 1 + x + x^2/2 + x^3/3 + x^4/12 }
  |>,
  <|
    "Label"->"Nonlinear exact (f=0,g=1)",
    "Eqns" -> { f'[x] == f[x] g[x], g'[x] == -f[x] },
    "ICs"  -> { f[0]==0, g[0]==1 },
    "Vars" -> {f, g}, "x0"->0, "Order"->6,
    "Expected" -> { 0, 1 }
  |>
};

(* ----------------------------- *)
(* Runner *)
(* ----------------------------- *)

printHeader["SYMMATH CROSS-CHECK OF MAXIMA TEST EXPECTATIONS (Mathematica)"];

(* Scalar ODEs *)
printSection["Scalar ODE Cases"];
scalarFailures = Reap[
  Do[
    Module[{lab, eqns, ics, var, x0, ord, exp, got, ok, arcsinHint},
      lab = c["Label"]; eqns = c["Eqns"]; ics = c["ICs"];
      var = c["Var"]; x0 = c["x0"]; ord = c["Order"]; exp = c["Expected"];
      arcsinHint = TrueQ[c["ArcsinHint"]];

      got = scalarSeriesSolve[eqns, ics, var, x, x0, ord];
      If[got === $Failed,
        printCase[lab, False, "FAILED to compute", exp, "ArcsinHint"->arcsinHint];
        Sow[lab];,
        ok = compare[got, exp];
        printCase[lab, ok, got, exp, "ArcsinHint"->arcsinHint];
        If[!ok, Sow[lab]];
      ];
    ],
    {c, scalarCases}
  ]
][[2, 1]] /. Null -> {};

(* 2×2 Systems *)
printSection["2×2 System Cases"];
systemFailures = Reap[
  Do[
    Module[{lab, eqns, ics, vars, x0, ord, expF, expG, gotList, gotF, gotG, okF, okG, hint},
      lab = s["Label"]; eqns = s["Eqns"]; ics = s["ICs"];
      vars = s["Vars"]; x0 = s["x0"]; ord = s["Order"]; hint = TrueQ[s["ExpSystemHint"]];
      {expF, expG} = s["Expected"];

      gotList = systemSeriesSolve2[eqns, ics, vars, x, x0, ord];
      If[gotList === $Failed || !MatchQ[gotList, {_, _}],
        printCase[lab, False, "FAILED to compute", {expF, expG}, "ExpSystemHint"->hint];
        Sow[lab];,
        {gotF, gotG} = gotList;
        okF = compare[gotF, expF];
        okG = compare[gotG, expG];
        If[okF && okG,
          printCase[lab, True, gotF, expF],
          (* When failing, print both f and g info *)
          Print["[FAIL] ", lab];
          Print["  (f) Got     : ", fmt[gotF]];
          Print["      Expected: ", fmt[expF]];
          Print["      Diff    : ", fmt[Expand[gotF - expF]]];
          Print["  (g) Got     : ", fmt[gotG]];
          Print["      Expected: ", fmt[expG]];
          Print["      Diff    : ", fmt[Expand[gotG - expG]]];
          If[hint,
            Print["  SUGGESTION: For f' = f + g, g' = f + g with f(0) = 1, g(0) = 0,"];
            Print["              f = (1 + E^(2 x))/2, g = (E^(2 x) - 1)/2 so x^5 coefficient is 2/15 for both,"];
            Print["              not 1/15. Update the expected series."];
          ];
          Sow[lab];
        ];
      ];
    ],
    {s, systemCases}
  ]
][[2, 1]] /. Null -> {};

(* Summary *)
Print["\n", StringRepeat["-", 72]];
If[scalarFailures === {} && systemFailures === {},
  Print["RESULT: ALL EXPECTATIONS MATCH MATHEMATICA (nice!)"],
  Print["RESULT: SOME TEST ASSUMPTIONS DISAGREE WITH MATHEMATICA"];
  If[scalarFailures =!= {},
    Print["Offending scalar cases:"]; Scan[Print["  - ", #] &, scalarFailures];
  ];
  If[systemFailures =!= {},
    Print["Offending system cases:"]; Scan[Print["  - ", #] &, systemFailures];
  ];
  Print["\nRecommended fixes (same as SymPy cross-check):"];
  Print["  1) Phase 3 Systems → \"Exp system f/g\": change x^5 coefficient from 1/15 to 2/15."];
  Print["  2) Special Functions → \"Integration gives arcsin\": "];
  Print["     EITHER change RHS to 1/Sqrt[1 - x^2], OR keep RHS Sqrt[1 - x^2)"];
  Print["     and change the expected polynomial to x - x^3/6 - x^5/40; also rename the test."];
];

(* Helpful tip mirroring your Maxima performance check issue *)
Print["\nPerformance test tip (Maxima): avoid length(string(res)) since 'length' expects a list."];
Print["  Use something robust instead, e.g.:"];
Print["      length(args(expand(res))) > 10   /* number of summed terms */"];
Print["    or hipow(res, x) >= 15             /* highest power check */"];
