(* File:    graph.wl    *)
(* Author:  Daniel Le   *)
(* Date:    07/10/2026  *)

(* Compute graph diameters and convert graph data structures to the string *)

BeginPackage["Graph`"]

graphDiameter::usage = "graphDiameter[graph_Association] returns {diameter, timesAttained}. \
graph must have keys \"n\" (vertex count), \"edges\" (list of 0-indexed {a,b} pairs, \
matching C++'s Edge set), and \"isDirected\" (True/False). Mirrors Graph::diameter, \
including its non-standard while(wasChange)-driven all-pairs shortest path sweep.";

graphToString::usage = "graphToString[graph_Association] returns a string matching \
Graph::toString's exact output format, including the per-edge trailing newlines and \
comma-on-its-own-line separators.";

Begin["`Private`"]

graphDiameter[graph_Association] := Module[
  {n, edges, isDirected, A, timesAttained2, ret},

  n = graph["n"];
  edges = graph["edges"];        (* List of 0-indexed {a,b} edge pairs *)
  isDirected = graph["isDirected"];

  (* Initialize the distance matrix. *)
  A = ConstantArray[n, {n, n}];

  (* Assign distance 1 to every edge. *)
  Do[
    (
      A[[edge[[1]] + 1, edge[[2]] + 1]] = 1;
      If[!isDirected,
        A[[edge[[2]] + 1, edge[[1]] + 1]]
      ] = 1
    ),
    {edge, edges}
  ];

  (* Set every vertex's distance to itself to 0. *)
  Do[
    A[[i, i]] = 0,
    {i, 1, n}
  ];

  (* Compute all-pairs shortest-path distances using the
     Floyd-Warshall algorithm. *)
  Do[
    A[[i, j]] = Min[A[[i, j]], A[[i, k]] + A[[k, j]]],
    {k, 1, n}, {i, 1, n}, {j, 1, n}
  ];

  (* Find the graph diameter (maximum finite shortest-path distance)
     and count the number of ordered vertex pairs attaining it. *)
  timesAttained2 = 0;
  ret = 0;

  Do[
    (
      If[A[[i, j]] > ret,
        ret = A[[i, j]];
        timesAttained2 = 0;
      ];
      If[A[[i, j]] == ret,
        timesAttained2++;
      ]
    ),
    {i, 1, n}, {j, 1, n}
  ];

  {ret, timesAttained2}
];

graphToString[graph_Association] := Module[{n, edges, m, lines, body},
  n = graph["n"];
  edges = graph["edges"];
  m = Length[edges];

  lines = Table["(" <> ToString[edge[[1]]] <> "," <> ToString[edge[[2]]] <> ")", {edge, edges}];

  body = If[m > 0, StringRiffle[lines, "\n,\n"], ""];

  If[m > 0,
    "(" <> ToString[n] <> ",{\n" <> body <> "\n}\n",
    "(" <> ToString[n] <> ",{\n}\n"
  ]
];

End[]
EndPackage[]
