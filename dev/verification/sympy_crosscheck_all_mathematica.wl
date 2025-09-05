(* ::Package:: *)

(*  sympy_crosscheck_all_mathematica.wl  (revised)
    ----------------------------------------------
    Run:
      "C:\Program Files\Wolfram Research\Wolfram\14.3\wolframscript.exe" -file sympy_crosscheck_all_mathematica.wl
*)

ClearAll["Global`*"];
$Assumptions = True;

fmt[expr_] := ToString[expr // Expand // Simplify, InputForm];

printHeader[title_] := (Print[StringRepeat["=",72]]; Print[title]; Print[StringRepeat["=",72]]);
printSection[title_] := Print["\n--- ", title, " ---"];

Options[printCase] = {"ArcsinHint"->False, "ExpSystemHint"->False};
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

compare[got_, expected_] := TrueQ[Simplify[Expand[got - expected] === 0]];

(* Extract initial value y[x0] == c from ICs; return $Failed if missing *)
initialValue[ics_List, y_Symbol, x0_] := Module[{hit},
  hit = Cases[ics, HoldPattern[y[x0] == v_] :> v, 1];
  If[hit === {}, $Failed, First@hit]
];

(* Scalar series with robust handling and ord==0 shortcut *)
scalarSeriesSolve[eqns_List, ics_List, y_, x_, x0_, ord_Integer] := Module[{ser, y0},
  If[ord == 0,
    y0 = initialValue[ics, y, x0];
    If[y0 =!= $Failed, Return[y0]];
  ];
  ser = Quiet@Check[
    Normal@AsymptoticDSolveValue[Join[eqns, ics], y[x], {x, x0, ord}],
    $Failed
  ];
  If[ser === $Failed,
    ser = Quiet@Check[
      With[{sol = DSolveValue[Join[eqns, ics], y[x], x]},
        Normal@Series[sol, {x, x0, ord}]
      ],
      $Failed
    ];
  ];
  ser
];

(* Manual series recursion for 2×2 systems f' = F(x,f,g), g' = G(x,f,g) *)
systemSeriesRecursion[Fexpr_, Gexpr_, f0_, g0_, x_, x0_, ord_Integer] := Module[
  {t, a, b, k, Ft, Gt, cF, cG, polyF, polyG},
  t = Unique["t"];
  a = Table[0, {ord+1}]; b = Table[0, {ord+1}];
  a[[1]] = f0; b[[1]] = g0;
  Do[
    Ft = Expand[
      Fexpr /. {x -> x0 + t,
                f[x] -> Sum[a[[i+1]] t^i, {i,0,k}],
                g[x] -> Sum[b[[i+1]] t^i, {i,0,k}]}
    ];
    Gt = Expand[
      Gexpr /. {x -> x0 + t,
                f[x] -> Sum[a[[i+1]] t^i, {i,0,k}],
                g[x] -> Sum[b[[i+1]] t^i, {i,0,k}]}
    ];
    cF = SeriesCoefficient[Normal@Series[Ft, {t,0,k}], {t,0,k}];
    cG = SeriesCoefficient[Normal@Series[Gt, {t,0,k}], {t,0,k}];
    a[[k+2]] = Simplify[cF/(k+1)];
    b[[k+2]] = Simplify[cG/(k+1)];
  , {k,0,ord-1}];
  polyF = Sum[a[[i+1]] (x - x0)^i, {i,0,ord}];
  polyG = Sum[b[[i+1]] (x - x0)^i, {i,0,ord}];
  {Expand@polyF, Expand@polyG}
];

systemSeriesSolve2[eqns_List, ics_List, vars_List, x_, x0_, ord_Integer] := Module[{Fexpr, Gexpr, f0, g0},
  If[vars =!= {f, g}, Return[$Failed]]; (* this script only supports {f,g} *)
  Fexpr = FirstCase[eqns, HoldPattern[f'[x] == rhs_] :> rhs, $Failed];
  Gexpr = FirstCase[eqns, HoldPattern[g'[x] == rhs_] :> rhs, $Failed];
  f0    = FirstCase[ics,  HoldPattern[f[x0] == c_]  :> c,   $Failed];
  g0    = FirstCase[ics,  HoldPattern[g[x0] == c_]  :> c,   $Failed];
  If[Fexpr === $Failed || Gexpr === $Failed || f0 === $Failed || g0 === $Failed, Return[$Failed]];
  systemSeriesRecursion[Fexpr, Gexpr, f0, g0, x, x0, ord]
];


xSym = x;

(* ----------------------------- *)
(* Build test cases (Scalar ODE) *)
(* ----------------------------- *)
scalarCases = {
  <|"Label"->"Linear homogeneous y'=y", "Eqns"->{y'[x]==y[x]}, "ICs"->{y[0]==1}, "Var"->y, "x0"->0, "Order"->6,
    "Expected"->1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 + x^6/720|>,
  <|"Label"->"Linear y'=2y", "Eqns"->{y'[x]==2 y[x]}, "ICs"->{y[0]==3}, "Var"->y, "x0"->0, "Order"->4,
    "Expected"->3 + 6 x + 6 x^2 + 4 x^3 + 2 x^4|>,
  <|"Label"->"Polynomial source y'=x^2", "Eqns"->{y'[x]==x^2}, "ICs"->{y[0]==0}, "Var"->y, "x0"->0, "Order"->5,
    "Expected"->x^3/3|>,
  <|"Label"->"Mixed polynomial y'=x+x^2", "Eqns"->{y'[x]==x + x^2}, "ICs"->{y[0]==1}, "Var"->y, "x0"->0, "Order"->4,
    "Expected"->1 + x^2/2 + x^3/3|>,
  <|"Label"->"Separable y'=y^2", "Eqns"->{y'[x]==y[x]^2}, "ICs"->{y[0]==1}, "Var"->y, "x0"->0, "Order"->5,
    "Expected"->1 + x + x^2 + x^3 + x^4 + x^5|>,
  <|"Label"->"Bernoulli y'=xy", "Eqns"->{y'[x]==x y[x]}, "ICs"->{y[0]==1}, "Var"->y, "x0"->0, "Order"->6,
    "Expected"->1 + x^2/2 + x^4/8 + x^6/48|>,
  <|"Label"->"Linear inhomo y'=y+sin(x)", "Eqns"->{y'[x]==y[x] + Sin[x]}, "ICs"->{y[0]==0}, "Var"->y, "x0"->0, "Order"->5,
    "Expected"->x^2/2 + x^3/6|>,
  <|"Label"->"Riccati y'=1+y^2", "Eqns"->{y'[x]==1 + y[x]^2}, "ICs"->{y[0]==0}, "Var"->y, "x0"->0, "Order"->5,
    "Expected"->x + x^3/3 + 2 x^5/15|>,
  With[{t=(x-1)},
    <|"Label"->"Nonzero point y'=x^2 at x=1", "Eqns"->{y'[x]==x^2}, "ICs"->{y[1]==2}, "Var"->y, "x0"->1, "Order"->4,
      "Expected"->2 + t + t^2 + t^3/3|>
  ],
  <|"Label"->"Exponential at x=2", "Eqns"->{y'[x]==y[x]}, "ICs"->{y[2]==1}, "Var"->y, "x0"->2, "Order"->4,
    "Expected"->1 + (x-2) + (x-2)^2/2 + (x-2)^3/6 + (x-2)^4/24|>,
  <|"Label"->"Linear at x=-1", "Eqns"->{y'[x]==2 x}, "ICs"->{y[-1]==0}, "Var"->y, "x0"->-1, "Order"->3,
    "Expected"->-2 (x+1) + (x+1)^2|>,
  <|"Label"->"Order 0", "Eqns"->{y'[x]==x + y[x]^2}, "ICs"->{y[0]==5}, "Var"->y, "x0"->0, "Order"->0,
    "Expected"->5|>,
  <|"Label"->"Order 1", "Eqns"->{y'[x]==3 x}, "ICs"->{y[0]==2}, "Var"->y, "x0"->0, "Order"->1,
    "Expected"->2|>,
  <|"Label"->"Constant RHS", "Eqns"->{y'[x]==5}, "ICs"->{y[0]==1}, "Var"->y, "x0"->0, "Order"->3,
    "Expected"->1 + 5 x|>,
  <|"Label"->"Complex coefficient", "Eqns"->{y'[x]==I y[x]}, "ICs"->{y[0]==1}, "Var"->y, "x0"->0, "Order"->3,
    "Expected"->1 + I x + (I^2) x^2/2 + (I^3) x^3/6|>,
  <|"Label"->"Hyperbolic cosine (y''=y)", "Eqns"->{y''[x]==y[x]}, "ICs"->{y[0]==1,y'[0]==0}, "Var"->y, "x0"->0, "Order"->6,
    "Expected"->1 + x^2/2 + x^4/24 + x^6/720|>,
  <|"Label"->"Trigonometric cosine (y''=-y)", "Eqns"->{y''[x]==-y[x]}, "ICs"->{y[0]==1,y'[0]==0}, "Var"->y, "x0"->0, "Order"->6,
    "Expected"->1 - x^2/2 + x^4/24 - x^6/720|>,
  <|"Label"->"Trigonometric sine (y''=-y, sine IC)", "Eqns"->{y''[x]==-y[x]}, "ICs"->{y[0]==0,y'[0]==1}, "Var"->y, "x0"->0, "Order"->5,
    "Expected"->x - x^3/6 + x^5/120|>,
  <|"Label"->"SHO frequency 2 (y''=-4y)", "Eqns"->{y''[x]==-4 y[x]}, "ICs"->{y[0]==1,y'[0]==0}, "Var"->y, "x0"->0, "Order"->4,
    "Expected"->1 - 2 x^2 + (2/3) x^4|>,
  <|"Label"->"Third order polynomial (y'''=0)", "Eqns"->{y'''[x]==0}, "ICs"->{y[0]==1,y'[0]==2,y''[0]==3}, "Var"->y, "x0"->0, "Order"->3,
    "Expected"->1 + 2 x + (3/2) x^2|>,
  <|"Label"->"Fourth order y^(4)=y", "Eqns"->{y''''[x]==y[x]}, "ICs"->{y[0]==1,y'[0]==0,y''[0]==0,y'''[0]==0}, "Var"->y, "x0"->0, "Order"->8,
    "Expected"->1 + x^4/24 + x^8/40320|>,
  <|"Label"->"Third order coupling (y'''=x+y)", "Eqns"->{y'''[x]==x + y[x]}, "ICs"->{y[0]==1,y'[0]==0,y''[0]==0}, "Var"->y, "x0"->0, "Order"->5,
    "Expected"->1 + x^3/6 + x^4/24|>,
  <|"Label"->"Second order at x=1 (y''=y)", "Eqns"->{y''[x]==y[x]}, "ICs"->{y[1]==1,y'[1]==1}, "Var"->y, "x0"->1, "Order"->4,
    "Expected"->1 + (x-1) + (x-1)^2/2 + (x-1)^3/6 + (x-1)^4/24|>,
  <|"Label"->"Constant third derivative at x=2 (y'''=6)", "Eqns"->{y'''[x]==6}, "ICs"->{y[2]==0,y'[2]==0,y''[2]==0}, "Var"->y, "x0"->2, "Order"->4,
    "Expected"->(x-2)^3|>,
  <|"Label"->"Integration of cosine", "Eqns"->{y'[x]==Cos[x]}, "ICs"->{y[0]==0}, "Var"->y, "x0"->0, "Order"->7,
    "Expected"->x - x^3/6 + x^5/120 - x^7/5040|>,
  <|"Label"->"Integration gives arctan", "Eqns"->{y'[x]==1/(1+x^2)}, "ICs"->{y[0]==0}, "Var"->y, "x0"->0, "Order"->9,
    "Expected"->x - x^3/3 + x^5/5 - x^7/7 + x^9/9|>,
  <|"Label"->"Integration gives arcsin (as written uses sqrt)", "Eqns"->{y'[x]==Sqrt[1 - x^2]}, "ICs"->{y[0]==0},
    "Var"->y, "x0"->0, "Order"->5, "Expected"->x + x^3/6 + 3 x^5/40, "ArcsinHint"->True|>
};

(* ----------------------------- *)
(* 2×2 systems *)
(* ----------------------------- *)
systemCases = {
  <|"Label"->"Simple oscillator f/g", "Eqns"->{f'[x]==g[x], g'[x]==-f[x]}, "ICs"->{f[0]==0, g[0]==1},
    "Vars"->{f,g}, "x0"->0, "Order"->5, "Expected"->{x - x^3/6 + x^5/120, 1 - x^2/2 + x^4/24}|>,
  <|"Label"->"Exp system f/g", "Eqns"->{f'[x]==f[x]+g[x], g'[x]==f[x]+g[x]}, "ICs"->{f[0]==1, g[0]==0},
    "Vars"->{f,g}, "x0"->0, "Order"->5,
    "Expected"->{1 + x + x^2 + (2 x^3)/3 + (x^4)/3 + (x^5)/15,
                 x + x^2 + (2 x^3)/3 + (x^4)/3 + (x^5)/15}, "ExpSystemHint"->True|>,
  <|"Label"->"Mixed poly f/g", "Eqns"->{f'[x]==x, g'[x]==f[x]+g[x]}, "ICs"->{f[0]==0, g[0]==1},
    "Vars"->{f,g}, "x0"->0, "Order"->4, "Expected"->{x^2/2, 1 + x + x^2/2 + x^3/3 + x^4/12}|>,
  <|"Label"->"Nonlinear exact (f=0,g=1)", "Eqns"->{f'[x]==f[x] g[x], g'[x]==-f[x]}, "ICs"->{f[0]==0, g[0]==1},
    "Vars"->{f,g}, "x0"->0, "Order"->6, "Expected"->{0,1}|>
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
  If[scalarFailures =!= {}, Print["Offending scalar cases:"]; Scan[Print["  - ", #] &, scalarFailures]];
  If[systemFailures =!= {}, Print["Offending system cases:"]; Scan[Print["  - ", #] &, systemFailures]];
  Print["\nRecommended fixes (same as SymPy cross-check):"];
  Print["  1) Phase 3 Systems → \"Exp system f/g\": change x^5 coefficient from 1/15 to 2/15."];
  Print["  2) Special Functions → \"Integration gives arcsin\": "];
  Print["     EITHER change RHS to 1/Sqrt[1 - x^2], OR keep RHS Sqrt[1 - x^2)"];
  Print["     and change the expected polynomial to x - x^3/6 - x^5/40; also rename the test."];
];

(* Maxima performance hint *)
Print["\nPerformance test tip (Maxima): avoid length(string(res)) since 'length' expects a list."];
Print["  Use something robust instead, e.g.:"];
Print["      length(args(expand(res))) > 10   /* number of summed terms */"];
Print["    or hipow(res, x) >= 15             /* highest power check */"];
