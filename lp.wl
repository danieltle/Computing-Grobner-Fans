BeginPackage["LpSolver`"]

(* Public Interface Functions *)
lpSetSolver::usage = "lpSetSolver[name_String] selects the current active LP solver.";
isFacet::usage = "isFacet[g_, i_]";
interiorPoint::usage = "interiorPoint[g_, result_, strictlyPositive_, equalitySet_]";
hasInteriorPoint::usage = "hasInteriorPoint[g_, strictlyPositive_, equalitySet_]";
shootRay::usage = "shootRay[g_]";
positiveVectorInKernel::usage = "positiveVectorInKernel[g_, result_]";
rankOfMatrix::usage = "rankOfMatrix[g_]";
extremeRaysInequalityIndices::usage = "extremeRaysInequalityIndices[inequalityList_]";
removeRedundantRows::usage = "removeRedundantRows[inequalities_, equalities_, removeInequalityRedundancies_]";
relativeInteriorPoint::usage = "relativeInteriorPoint[n_, g_, equalitySet_]";
dual::usage = "dual[n_, inequalities_, equations_, dualInequalities_, dualEquations_]";
hasHomogeneousSolution::usage = "hasHomogeneousSolution[n_, inequalities_, equations_]";
isInNonNegativeSpan::usage = "isInNonNegativeSpan[v_, rays_, linealitySpace_]";

(* Simulation of the abstract base class and registration registry *)
CreateLpSolver::usage = "CreateLpSolver[name, methodsAssoc] registers a new solver instance.";
PrintSolverList::usage = "PrintSolverList[] prints all registered LP solvers.";

Begin["`Private`"]

(* Global Mutable States *)
$SolverList = <||>;
$DefaultSolver = None;
$Initialized = False;


(* Helper: primitive integer row reduction *)
primitiveRow[row_List] := Module[{cleared, g},
  (* Clear denominators first so GCD is taken over integers *)
  cleared = row * LCM @@ (Denominator /@ row);
  g = GCD @@ (Numerator /@ cleared);
  If[g == 0, cleared, cleared / g]
];

primitiveMatrix[m_List] := primitiveRow /@ m;

(* Constructor / Registry for Solvers *)
CreateLpSolver[name_String, methods_Association] := (
  $SolverList[name] = <|
    "Name" -> name,
    "interiorPoint" -> (MessageDialog["interiorPoint method not supported in " <> name]; Abort[] &),
    "hasInteriorPoint" -> (MessageDialog["hasInteriorPoint method not supported in " <> name]; Abort[] &),
    "shoot" -> (MessageDialog["shoot method not supported in " <> name]; Abort[] &),
    "positiveVectorInKernel" -> (MessageDialog["positiveVectorInKernel method not supported in " <> name]; Abort[] &),
    "rankOfMatrix" -> (MessageDialog["rankOfMatrix method not supported in " <> name]; Abort[] &),
    "extremeRaysInequalityIndices" -> (MessageDialog["extremeRaysInequalityIndices not supported in " <> name]; Abort[] &),
    "removeRedundantRows" -> (MessageDialog["removeRedundantRows method not supported in " <> name]; Abort[] &),
    "relativeInteriorPoint" -> (MessageDialog["relativeInteriorPoint method not supported in " <> name]; Abort[] &),
    "dual" -> (MessageDialog["dual method not supported in " <> name]; Abort[] &),
    "hasHomogeneousSolution" -> (MessageDialog["hasHomogeneousSolution method not supported in " <> name]; Abort[] &),
    (* FIX: isFacet no longer silently defaults to True. C++'s base   *)
    (* class has no concrete isFacet implementation (every subclass   *)
    (* must override it), so the fallback here now matches every      *)
    (* other method's "not supported" abort behavior instead of       *)
    (* fabricating a result.                                          *)
    "isFacet" -> (MessageDialog["isFacet method not supported in " <> name]; Abort[] &)
  |>;
  (* Merge user-provided overridden concrete implementations *)
  $SolverList[name] = Merge[{$SolverList[name], methods}, Last];
  $SolverList[name]
);

PrintSolverList[] := (
  Print["List of linked LP solvers:"];
  Scan[Print[" ", #["Name"]] &, Values[$SolverList]]
);

(* Initialization Logic *)
LpInit[] := If[!$Initialized, lpSetSolver[""]];

lpSetSolver[name_String] := Module[{selected, soplexCddGmp, soplex, huber, cdd, cddgmp},
  soplexCddGmp = Lookup[$SolverList, "SoPlexCddGmp", None];
  soplex = Lookup[$SolverList, "SoPlex", None];
  huber = Lookup[$SolverList, "Huber's", None];
  cdd = Lookup[$SolverList, "cdd", None];
  cddgmp = Lookup[$SolverList, "cddgmp", None];
  selected = Lookup[$SolverList, name, None];

  (* Cascading fallback initialization logic identical to C++ *)
  $DefaultSolver = huber;
  If[soplex =!= None, $DefaultSolver = soplex];
  If[cdd =!= None, $DefaultSolver = cdd];
  If[cddgmp =!= None, $DefaultSolver = cddgmp];
  If[soplexCddGmp =!= None, $DefaultSolver = soplexCddGmp];
  If[selected =!= None, $DefaultSolver = selected];

  $Initialized = True;

  If[$DefaultSolver === None,
    Print["Error: No default LP solver could be asserted."]; Abort[],
    Print["LP algorithm being used: \"", $DefaultSolver["Name"], "\"."]
  ];
  selected
];

(* Interface Wrapper Implementations *)

isFacet[g_, i_] := (LpInit[]; $DefaultSolver["isFacet"][g, i]);

interiorPoint[g_, result_, strictlyPositive_, equalitySet_] := (LpInit[]; $DefaultSolver["interiorPoint"][g, result, strictlyPositive, equalitySet]);

hasInteriorPoint[g_, strictlyPositive_, equalitySet_] := (LpInit[]; $DefaultSolver["hasInteriorPoint"][g, strictlyPositive, equalitySet]);

shootRay[g_] := (
  LpInit[];
  If[Length[g] == 0, Return[None]];
  $DefaultSolver["shoot"][g]
);

positiveVectorInKernel[g_, result_] := (LpInit[]; $DefaultSolver["positiveVectorInKernel"][g, result]);

rankOfMatrix[g_] := (LpInit[]; $DefaultSolver["rankOfMatrix"][g]);

extremeRaysInequalityIndices[inequalityList_] := Module[{m, ret},
  (* Simplicial cone optimization check *)
  If[rankOfMatrix[inequalityList] == Length[inequalityList],
    m = Length[inequalityList];
    ret = Table[
      Join[Range[1, i - 1], Range[i + 1, m]],
      {i, 1, m}
    ];
    Return[ret]
  ];
  LpInit[];
  $DefaultSolver["extremeRaysInequalityIndices"][inequalityList]
];

removeRedundantRows[inequalities_, equalities_, removeInequalityRedundancies_] := (
  LpInit[];
  $DefaultSolver["removeRedundantRows"][inequalities, equalities, removeInequalityRedundancies]
);

relativeInteriorPoint[n_, g_, equalitySet_] := (LpInit[]; $DefaultSolver["relativeInteriorPoint"][n, g, equalitySet]);

dual[n_, inequalities_, equations_, dualInequalities_, dualEquations_] := (
  LpInit[];
  $DefaultSolver["dual"][n, inequalities, equations, dualInequalities, dualEquations]
);

hasHomogeneousSolution[n_Integer, inequalities_List, equations_List] := (
  LpInit[];
  (* Input integrity validation *)
  Scan[
    If[Length[#] != n,
      Print["Inequality length does not match. n=", n, " vector=", #];
      Abort[]
    ] &,
    inequalities
  ];
  Scan[
    If[Length[#] != n,
      Print["Equation length does not match. n=", n, " vector=", #];
      Abort[]
    ] &,
    equations
  ];

  $DefaultSolver["hasHomogeneousSolution"][n, inequalities, equations]
);

(* Linear Algebra Logic Integration *)
isInNonNegativeSpan[v_List, rays_List, linealitySpace_List] := Module[
  {n, numRays, numLin, A1, temp, A2, PrimitiveMatrixA1, PrimitiveMatrixA2},

  n = Length[v];
  numRays = Length[rays];
  numLin = Length[linealitySpace];

  (* Constructing A1 Matrix block components *)
  A1 = Join[
    ConstantArray[0, {numRays, 1}],
    IdentityMatrix[numRays],
    ConstantArray[0, {numRays, numLin}],
    2
  ];

  (* Constructing A2. C++ builds A2 by horizontally transposed pieces *)
  temp = {-v}; (* 1 x n matrix *)
  A2 = Transpose[Join[temp, rays, linealitySpace]]; (* n x (1+numRays+numLin) *)
  
  PrimitiveMatrixA1 = primitiveMatrix[A1];
  PrimitiveMatrixA2 = primitiveMatrix[A2];

  hasHomogeneousSolution[
    1 + numRays + numLin,
    PrimitiveMatrixA1,
    PrimitiveMatrixA2
  ]
];

End[]
EndPackage[]
