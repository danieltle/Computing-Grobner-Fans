(* File:    termorder.wl   *)
(* Author:  Daniel Le      *)
(* Date:    07/11/2026     *)

(* Defines term orders *)


BeginPackage["TermOrder`"]

TermOrderRowDot::usage = "TermOrderRowDot[to_, row_, v_]";
TermOrderCompare::usage = "TermOrderCompare[to_, a_, b_, scaleA_:1, scaleB_:1, perturbationDegree_:-1]";
TermOrderPrint::usage = "TermOrderPrint[to_] mirrors TermOrder::print.";
TermOrderPrintMatrix::usage = "TermOrderPrintMatrix[to_, dim_] mirrors TermOrder::printMatrix.";

LexicographicTermOrder::usage = "LexicographicTermOrder[largest_:0]";
LexicographicInvertedTermOrder::usage = "LexicographicInvertedTermOrder[]";
ReverseLexicographicTermOrder::usage = "ReverseLexicographicTermOrder[largest_:0]";
ReverseLexicographicInvertedTermOrder::usage = "ReverseLexicographicInvertedTermOrder[]";
StandardGradedLexicographicTermOrder::usage = "StandardGradedLexicographicTermOrder[largest_:0] ";
WeightTermOrder::usage = "WeightTermOrder[weight_]";
WeightReverseLexicographicTermOrder::usage = "WeightReverseLexicographicTermOrder[weight_]";
MatrixTermOrder::usage = "MatrixTermOrder[weights_] where weights is a list of IntegerVectors (rows).";
TotalDegreeTieBrokenTermOrder::usage = "TotalDegreeTieBrokenTermOrder[tieBreaker_] where tieBreaker is another term order Association.";

Begin["`Private`"]

TermOrderRowDot[to_Association, row_Integer, v_List] := to["RowDotImpl"][row, v];

TermOrderCompare[to_Association, a_List, b_List, scaleA_: 1, scaleB_: 1, perturbationDegree_: -1] :=
  to["CompareImpl"][a, b, scaleA, scaleB, perturbationDegree];

TermOrderPrint[to_Association] := to["PrintImpl"][];

(* TermOrder::printMatrix - uses standardVector(dim,j), the j-th standard basis vector of length dim *)
standardVector[dim_Integer, j_Integer] := ReplacePart[ConstantArray[0, dim], j + 1 -> 1];

TermOrderPrintMatrix[to_Association, dim_Integer] := Module[{l},
  l = Table[
    Table[TermOrderRowDot[to, i, standardVector[dim, j]], {j, 0, dim - 1}],
    {i, 0, dim - 1}
  ];
  Print[l];
  Print[""];
];

dot[v_List, w_List] := v.w;
dotLong[v_List, w_List] := v.w;

(* LexicographicTermOrder *)

LexicographicTermOrder[largest_Integer: 0] := <|
  "Name" -> "LexicographicTermOrder",
  "Largest" -> largest,

  "RowDotImpl" -> Function[{row, v},
    v[[Mod[row + largest, Length[v]] + 1]]
  ],

  "CompareImpl" -> Function[{a, b, scaleA, scaleB, perturbationDegree},
    Module[{n, nLoop, idx},
      If[Length[a] != Length[b],
        Print["Lexicographic term order compare failed on the following vectors:"];
        Print[a];
        Print[b];
        Abort[]
      ];
      n = Length[a];
      nLoop = If[perturbationDegree >= 0, perturbationDegree, n];
      Do[
        idx = Mod[i + largest, n] + 1;
        If[scaleA*a[[idx]] < scaleB*b[[idx]], Return[True]];
        If[scaleA*a[[idx]] > scaleB*b[[idx]], Return[False]];
        ,
        {i, 0, nLoop - 1}
      ];
      False
    ]
  ],

  "PrintImpl" -> (Print["LexicographicTermOrder"]; Print[""]; &)
|>;

(* LexicographicInvertedTermOrder *)

LexicographicInvertedTermOrder[] := <|
  "Name" -> "LexicographicInvertedTermOrder",

  "RowDotImpl" -> Function[{row, v},
    v[[Length[v] - row]]  (* C++: v[v.size()-row-1], 0-indexed -> +1 for Wolfram *)
  ],

  "CompareImpl" -> Function[{a, b, scaleA, scaleB, perturbationDegree},
    Module[{n},
      n = Length[a];
      Do[
        If[scaleA*a[[i + 1]] < scaleB*b[[i + 1]], Return[True]];
        If[scaleA*a[[i + 1]] > scaleB*b[[i + 1]], Return[False]];
        ,
        {i, n - 1, 0, -1}
      ];
      False
    ]
  ],

  "PrintImpl" -> (Print["LexicographicInvertedTermOrder"]; Print[""]; &)
|>;

(* ReverseLexicographicTermOrder *)

ReverseLexicographicTermOrder[largest_Integer: 0] := <|
  "Name" -> "ReverseLexicographicTermOrder",
  "Largest" -> largest,

  (* index(row,a) = (-row + a.size() + largest - 1) mod a.size(), 0-indexed *)
  "IndexImpl" -> Function[{row, a}, Mod[-row + Length[a] + largest - 1, Length[a]]],

  "RowDotImpl" -> Function[{row, v},
    -v[[Mod[-row + Length[v] + largest - 1, Length[v]] + 1]]
  ],

  "CompareImpl" -> Function[{a, b, scaleA, scaleB, perturbationDegree},
    Module[{n, nLoop, idx, A, B},

      If[Length[a] != Length[b], Abort[]];
      n = Length[a];
      nLoop = If[perturbationDegree >= 0, perturbationDegree, n];
      Do[
        idx = Mod[-i + n + largest - 1, n] + 1;
        A = scaleA*a[[idx]];
        B = scaleB*b[[idx]];
        If[A > B, Return[True]];
        If[A < B, Return[False]];
        ,
        {i, 0, nLoop - 1}
      ];
      False
    ]
  ],

  "PrintImpl" -> (Print["ReverseLexicographicTermOrder"]; Print[""]; &)
|>;

(* ReverseLexicographicInvertedTermOrder *)

ReverseLexicographicInvertedTermOrder[] := <|
  "Name" -> "ReverseLexicographicInvertedTermOrder",

  "RowDotImpl" -> Function[{row, v}, -v[[row + 1]]],

  "CompareImpl" -> Function[{a, b, scaleA, scaleB, perturbationDegree},
    Module[{n, A, B},
      (* Bare assert, no message, matching C++. *)
      If[Length[a] != Length[b], Abort[]];
      n = Length[a];
      Do[
        A = scaleA*a[[i + 1]];
        B = scaleB*b[[i + 1]];
        If[A != B,
          If[A > B, Return[True]];
          If[A < B, Return[False]];
        ];
        ,
        {i, 0, n - 1}
      ];
      False
    ]
  ],

  "PrintImpl" -> (Print["ReverseLexicographicInvertedTermOrder"]; Print[""]; &)
|>;

(* StandardGradedLexicographicTermOrder *)

StandardGradedLexicographicTermOrder[largest_Integer: 0] := <|
  "Name" -> "StandardGradedLexicographicTermOrder",
  "Largest" -> largest,

  "RowDotImpl" -> Function[{row, v},
    If[row == 0,
      Total[v],
      v[[Mod[(row - 1) + largest, Length[v]] + 1]]
    ]
  ],

  "CompareImpl" -> Function[{a, b, scaleA, scaleB, perturbationDegree},
    Module[{n, difsum, nLoop, idx},
      (* Bare assert, no message, matching C++. *)
      If[Length[a] != Length[b], Abort[]];
      n = Length[a];

      If[perturbationDegree == 0, Return[False]];

      difsum = scaleA*Total[a] - scaleB*Total[b];
      If[difsum < 0, Return[True]];
      If[difsum > 0, Return[False]];

      If[perturbationDegree == 1, Return[False]];

      nLoop = If[perturbationDegree >= 0, perturbationDegree - 1, n];
      Do[
        idx = Mod[i + largest, n] + 1;
        If[scaleA*a[[idx]] < scaleB*b[[idx]], Return[True]];
        If[scaleA*a[[idx]] > scaleB*b[[idx]], Return[False]];
        ,
        {i, 0, nLoop - 1}
      ];
      False
    ]
  ],

  "PrintImpl" -> (Print["StandardGradedLexicographicTermOrder"]; Print[""]; &)
|>;

(* WeightTermOrder *)

WeightTermOrder[weight_List] := <|
  "Name" -> "WeightTermOrder",
  "Weight" -> weight,

  "RowDotImpl" -> Function[{row, v},
    If[row == 0,
      dot[v, weight],
      TermOrderRowDot[LexicographicTermOrder[], row - 1, v]  (* default largest=0, see usage note *)
    ]
  ],

  "CompareImpl" -> Function[{a, b, scaleA, scaleB, perturbationDegree},
    Module[{d},
      If[perturbationDegree == 0, Return[False]];
      d = scaleA*dotLong[a, weight] - scaleB*dotLong[b, weight];
      If[d < 0, Return[True]];
      If[d > 0, Return[False]];
      TermOrderCompare[LexicographicTermOrder[], a, b, scaleA, scaleB, perturbationDegree - 1]
    ]
  ],

  "PrintImpl" -> (Print["WeightTermOrder"]; Print[weight]; Print[""]; &)
|>;

(* WeightReverseLexicographicTermOrder *)

WeightReverseLexicographicTermOrder[weight_List] := <|
  "Name" -> "WeightReverseLexicographicTermOrder",
  "Weight" -> weight,

  "GetWeightImpl" -> Function[{}, weight],

  "RowDotImpl" -> Function[{row, v},
    If[row == 0,
      dot[v, weight],
      TermOrderRowDot[ReverseLexicographicTermOrder[], row - 1, v]  (* default largest=0 *)
    ]
  ],

  "CompareImpl" -> Function[{a, b, scaleA, scaleB, perturbationDegree},
    Module[{d},
      If[perturbationDegree == 0, Return[False]];
      d = scaleA*dotLong[a, weight] - scaleB*dotLong[b, weight];
      If[d < 0, Return[True]];
      If[d > 0, Return[False]];
      TermOrderCompare[ReverseLexicographicTermOrder[], a, b, scaleA, scaleB, perturbationDegree - 1]
    ]
  ],

  "PrintImpl" -> (Print["WeightReverseLexicographicTermOrder"]; Print[weight]; Print[""]; &)
|>;

(* MatrixTermOrder *)

MatrixTermOrder[weights_List] := <|
  "Name" -> "MatrixTermOrder",
  "Weights" -> weights,

  "RowDotImpl" -> Function[{row, v},
    Module[{nrows = Length[weights]},
      If[row < nrows,
        dot[v, weights[[row + 1]]],
        TermOrderRowDot[ReverseLexicographicInvertedTermOrder[], row - nrows, v]
      ]
    ]
  ],

  "CompareImpl" -> Function[{a, b, scaleA, scaleB, perturbationDegree},
    Module[{da, db, d, n},
    
      If[perturbationDegree == 0, Abort[]];

      n = Length[a];
      Do[
        da = weightRow.a;
        db = weightRow.b;
        d = scaleA*da - scaleB*db;
        If[d < 0, Return[True]];
        If[d > 0, Return[False]];
        ,
        {weightRow, weights}
      ];

      TermOrderCompare[ReverseLexicographicInvertedTermOrder[], a, b, scaleA, scaleB, perturbationDegree - 1]
    ]
  ],

  "PrintImpl" -> (Print["MatrixTermOrder"]; Print[weights]; Print[""]; &)
|>;

(* TotalDegreeTieBrokenTermOrder *)

TotalDegreeTieBrokenTermOrder[tieBreaker_Association] := <|
  "Name" -> "TotalDegreeTieBrokenTermOrder",
  "TieBreaker" -> tieBreaker,

  "RowDotImpl" -> Function[{row, v},
    If[row == 0,
      Total[v],
      TermOrderRowDot[tieBreaker, row - 1, v]
    ]
  ],

  "CompareImpl" -> Function[{a, b, scaleA, scaleB, perturbationDegree},
    Module[{d},
      d = scaleA*Total[a] - scaleB*Total[b];
      If[d < 0, Return[True]];
      If[d > 0, Return[False]];
      TermOrderCompare[tieBreaker, a, b, scaleA, scaleB, perturbationDegree - 1]
    ]
  ],

  "PrintImpl" -> (Print["TotalDegreeTieBrokenTermOrder"]; TermOrderPrint[tieBreaker]; Print[""]; &)
|>;

End[]
EndPackage[]
