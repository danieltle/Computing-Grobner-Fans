(*1. BUILD WEIGHT MATRIX FOR MONOMIAL ORDER*)(*-------------------------------------------------------------*)(*C\
onvert a weight vector to a square weight matrix for MonomialOrder*)
weightMatrix[w_List] := 
 Module[{n = Length[w]},(*First row:
  the weight vector*)(*Remaining rows:identity matrix for tie-
  breaking*)Join[{w}, IdentityMatrix[n]][]]

(*Alternative:with degree-aware tie-breaker*)
weightMatrixWithTieBreaker[w_List, 
  tieBreaker_String : "Lexicographic"] := 
 Module[{n = Length[w]}, 
  Switch[tieBreaker, "Lexicographic", Join[{w}, IdentityMatrix[n]], 
   "DegreeReverseLexicographic", 
   Join[{w}, Table[If[i + j == n + 1, 1, 0], {i, n}, {j, n}]], _, 
   Join[{w}, IdentityMatrix[n]]]]

(*-------------------------------------------------------------*)
(*2. EXTRACT INITIAL FORM WITH RESPECT TO A WEIGHT VECTOR*)
(*-------------------------------------------------------------*)

initialForm[expr_, vars_, w_] := 
 Module[{terms, weights, maxW}, terms = List @@ Expand[expr];
  If[terms === {}, Return[0]];
  weights = w . Exponent[#, vars] & /@ terms;
  maxW = Max[weights];
  Plus @@ Pick[terms, weights, maxW]]

(*-------------------------------------------------------------*)
(*3. FIND THE FIRST FACET CROSSED ALONG THE SEGMENT*)
(*-------------------------------------------------------------*)

findNextFacet[G_List, vars_List, wCurr_List, wTarget_List] := 
 Module[{minT = Infinity, bestF = Null, bestBeta = Null, terms, 
   leadTerm, leadExp, beta, num, den, t}, Do[terms = List @@ Expand[f];
   If[terms === {}, Continue[]];
   (*Leading term at current weight*)
   leadTerm = First@SortBy[terms, -wCurr . Exponent[#, vars] &];
   leadExp = Exponent[leadTerm, vars];
   (*Compare all other terms to find where the leading term changes*)
   Do[beta = Exponent[term, vars];
    If[beta != leadExp, num = wCurr . (beta - leadExp);
     den = (wTarget - wCurr) . (beta - leadExp);
     (*t=num/den;require 0<t<1*)If[den > 0 && num > 0, t = num/den;
      If[t < minT && t < 1, minT = t;
       bestF = f;
       bestBeta = beta]]], {term, terms}], {f, G}];
  {minT, bestF, bestBeta}]

(*-------------------------------------------------------------*)
(*4. PERFORM A LOCAL FLIP ACROSS A FACET (Algorithm 4.2)*)
(*-------------------------------------------------------------*)

performFlip[G_List, vars_List, wCurr_List, wCross_List, 
  wTarget_List] := 
 Module[{G1, wNext, G2, lifted, newGB, n, wmat}, n = Length[vars];
  (*Step 1:Compute initial forms at the crossing point*)
  G1 = initialForm[#, vars, wCross] & /@ G;
  G1 = Select[G1, # =!= 0 &];
  If[G1 === {}, Return[G]];
  (*Step 2:
  Move infinitesimally past the facet*)(*For exact arithmetic,
  we use a small rational perturbation*)
  wNext = wCross + (1/10000) (wTarget - wCurr);
  (*Build weight matrix for the new side*)wmat = weightMatrix[wNext];
  (*Compute the other reduced Gröbner basis of the initial ideal*)
  G2 = GroebnerBasis[G1, vars, MonomialOrder -> wmat];
  G2 = Select[G2, # =!= 0 &];
  If[G2 === {}, Return[G]];
  (*Step 3:Lift back to the original ideal I*)(*For each g in G2,
  subtract its normal form modulo G*)
  lifted = Table[g2 - PolynomialReduce[g2, G, vars][[1]] . G, {g2, G2}];
  (*Remove zero polynomials*)newGB = Select[lifted, Expand[#] =!= 0 &];
  If[newGB === {}, Return[G]];
  (*Return the new basis*)newGB]

(*-------------------------------------------------------------*)
(*5. MAIN GENERIC GRÖBNER WALK FUNCTION*)
(*-------------------------------------------------------------*)

Options[groebnerWalk] = {MaxIterations -> 1000, Verbose -> True, 
   Epsilon -> 1/10000};

groebnerWalk[G0_List, vars_List, w0_List, wTarget_List, 
  opts : OptionsPattern[]] := 
 Module[{G = G0, wCurr = w0, n = Length[vars], 
   maxIter = OptionValue[MaxIterations], 
   verbose = OptionValue[Verbose], eps = OptionValue[Epsilon], 
   iter = 0, res, t, wCross, wmat, newG}, 
  If[verbose, Print["Starting Gröbner walk"];
   Print["  Initial basis: ", G];
   Print["  Start weight: ", w0];
   Print["  Target weight: ", wTarget];
   Print[""];];
  While[iter < maxIter,(*Find the first facet crossing*)
   res = findNextFacet[G, vars, wCurr, wTarget];
   {t, f, beta} = res;
   (*If no crossing before t=1,we're done*)
   If[t >= 1 || t === Infinity, 
    If[verbose, Print["No more facets to cross. Reached target cone."]];
    Break[]];
   (*Weight vector at the crossing*)
   wCross = wCurr + t (wTarget - wCurr);
   If[verbose, Print["Iteration ", iter + 1];
    Print["  Crossing facet at t = ", N[t]];
    Print["  Crossing weight: ", wCross];];
   (*Perform the flip to the neighboring cone*)
   G = performFlip[G, vars, wCurr, wCross, wTarget];
   If[verbose, Print["  New basis size: ", Length[G]];];
   (*Update current weight to just past the facet*)
   wCurr = wCross + eps (wTarget - wCurr);
   iter++;];
  If[verbose, Print[""];
   Print["Walk completed after ", iter, " iterations"];
   Print["Final basis size: ", Length[G]];];
  (*Compute the reduced Gröbner basis at the target weight*)
  wmat = weightMatrix[wTarget];
  GroebnerBasis[G, vars, MonomialOrder -> wmat]]

(*-------------------------------------------------------------*)
(*6. TEST EXAMPLE FROM THE PAPER*)
(*-------------------------------------------------------------*)

testGroebnerWalk[] := 
 Module[{vars = {x, y, 
     z},(*Gröbner basis for lex order from Example 2.7*)
   G0 = {y^2 + x - x^3 y - x^4, z + y + x},(*Weight for lex order:(1,
   4,5)*)w0 = {1, 4, 5},(*Target:degree reverse lex (approximately (1,
   1,1))*)wTarget = {1, 1, 1}, result}, Print["=" . 60];
  Print["TEST: Computing Gröbner walk for Example 2.7"];
  Print["=" . 60];
  Print["Ideal generators: ", {y^2 + x - x^3 y - x^4, z + y + x}];
  Print["Starting weight (lex): ", w0];
  Print["Target weight: ", wTarget];
  Print[""];
  result = groebnerWalk[G0, vars, w0, wTarget, Verbose -> True];
  Print[""];
  Print["Final reduced Gröbner basis:"];
  Print[result];
  Print[""];
  Print["Verification: The basis should be a Gröbner basis"];
  Print["for the term order with weight ", wTarget];
  result]

(*Uncomment to run the test*)
(*testGroebnerWalk[]*)

(*-------------------------------------------------------------*)
(*7. ADDITIONAL UTILITY:COMPARE MONOMIALS WITH WEIGHT*)
(*-------------------------------------------------------------*)

monomialCompare[m1_List, m2_List, w_List] := 
 Module[{diff = m1 - m2}, w . diff]

(*Check if a polynomial is homogeneous with respect to a weight*)
isHomogeneous[expr_, vars_, w_] := 
 Module[{terms = List @@ Expand[expr], weights}, 
  If[terms === {}, True];
  weights = w . Exponent[#, vars] & /@ terms;
  Length[Union[weights]] == 1]

(*-------------------------------------------------------------*)
(*8. EXAMPLE:PRINCIPAL IDEAL FROM EXAMPLE 2.9*)
(*-------------------------------------------------------------*)

testPrincipalIdeal[] := 
 Module[{vars = {x, y}, f = x^4 + x^4 y - x^3 y + x^2 y^2 + y, G0, w0,
    wTarget, result}, Print["Principal ideal: <", f, ">"];
  (*Initial basis is just the generator*)G0 = {f};
  (*Start with lex order:x>y*)w0 = {1, 0};
  wTarget = {0, 1};
  result = 
   groebnerWalk[G0, vars, w0, wTarget, Verbose -> True, 
    MaxIterations -> 10];
  Print["Final basis for y > x: ", result];
  result]

(*Uncomment to run*)
(*testPrincipalIdeal[]*)

(*-------------------------------------------------------------*)
(*9. VISUALIZATION AID:SHOW THE WALK PATH*)
(*-------------------------------------------------------------*)

traceGroebnerWalk[G0_List, vars_List, w0_List, wTarget_List] := 
 Module[{G = G0, wCurr = w0, n = Length[vars], path = {}, iter = 0, 
   maxIter = 100, res, t, wCross, wmat}, AppendTo[path, {w0, Length[G]}];
  While[iter < maxIter, res = findNextFacet[G, vars, wCurr, wTarget];
   {t, f, beta} = res;
   If[t >= 1 || t === Infinity, Break[]];
   wCross = wCurr + t (wTarget - wCurr);
   G = performFlip[G, vars, wCurr, wCross, wTarget];
   wCurr = wCross + (1/10000) (wTarget - wCurr);
   iter++;
   AppendTo[path, {wCurr, Length[G]}];];
  Print["Walk path:"];
  Print["Step\tWeight\t\tBasis size"];
  Do[Print[step, "\t", Round[path[[step, 1]], 0.001], "\t", 
    path[[step, 2]]], {step, Length[path]}];
  path]

(*-------------------------------------------------------------*)
(*10. END OF SKELETON*)
(*-------------------------------------------------------------*)

Print["Gröbner walk skeleton loaded successfully."];
Print["Functions available:"];
Print["  groebnerWalk[G, vars, w0, wTarget] - Main walk function"];
Print["  testGroebnerWalk[] - Test with example from paper"];
Print["  traceGroebnerWalk[G, vars, w0, wTarget] - Show walk path"];
Print["  weightMatrix[w] - Convert weight vector to matrix"];