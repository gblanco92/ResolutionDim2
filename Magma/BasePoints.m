import "ProximityMatrix.m": ProximityMatrixImpl;

intrinsic BasePoints(I::RngMPol : Coefficients := false) -> []
{ Computes the weighted cluster of base points of a bivariate
  polynomial ideal I }
require Type(Representative(I)) eq RngMPolElt:
  "Ideal must be a polynomial ideal";
require Rank(Parent(Representative(I))) eq 2:
  "Ideal must be a bivariate polynomial ideal";
  // Generators in G & fixed part F.
  G := Generators(I); F := Gcd(G);
  //G := [g div F : g in G] cat (Evaluate(F, <0,0>) eq 0 select [F] else []);
  G := [g div F : g in G] cat [F];
  // Compute the Puiseux expansion & the prox. matrix for the product.
  S := NewtonPuiseuxAlgorithm(G: Polynomial := true);
  P := ProximityMatrixImpl([* <s[1], 1> : s in S *]: ExtraPoint := true,
    Coefficients := true);
  C := P[3]; EE := P[2]; n := Ncols(P[1]); P := P[1]; ZZ := IntegerRing();
  // Compute the multiplicities of each generator in G.
  E := [ZeroMatrix(ZZ, 1, n) : i in [1..#G]];
  for i in [1..#S] do
    for m in S[i][2] do E[m[2]] := E[m[2]] + m[1] * EE[i]; end for;
  end for;
  // Values for each generator in G.
  V := [e*Transpose(P^-1) : e in E];
  // Values for the points in the cluster of base points.
  v := ZeroMatrix(ZZ, 1, n);
  for i in [1..n] do v[1][i] := Min([vj[1][i] : vj in Prune(V)]); end for;
  // Multiplicities for the cluster of base points.
  e := v*Transpose(P);
  // Remove points not in the cluster of base points.
  I := [i : i in [1..n] | e[1][i] ne 0 or E[#E][1][i] ne 0];
  P := Submatrix(P, I, I);  v := Submatrix(v, [1], I); n := Ncols(P);
  V := [Submatrix(v, [1], I) : v in V]; E := [Submatrix(e, [1], I) : e in E];
  // ------------ Add NEW free points ------------------
  e := v*Transpose(P); inCluster := [i : i in [1..n] | e[1][i] ne 0];
  S := &+Submatrix(P, inCluster, inCluster)[1..#inCluster];
  lastFree := [i : i in [1..#inCluster] | S[i] eq 1];
  // For each last free point on a branch...
  for p in lastFree do
    // Values for each generator in the point p.
    Vp := [vi[1][p] : vi in Prune(V)];
    // Index of the generator achieving the minimum.
    g := Index(Vp, Min(Vp));
    // If there is a unique gen. achieving the min. and
    // the excess is positive add new points.
    uniqueGen := #[vp : vp in Vp | vp eq v[1][p]] eq 1;
    if uniqueGen and E[g][1][p] ne 0 then
      // Minimm of the values for all the generator except g.
      wp := Min(Exclude(Vp, v[1][p]));
      // Number of new free points.
      k := Ceiling((wp - v[1][p])/E[g][1][p]);
      // Expand the proximity matrix.
      P := InsertBlock(ScalarMatrix(n + k, 1), P, 1, 1);
      P[n + 1][p] := -1; n := Ncols(P);
      for i in [1..k - 1] do P[n - i + 1][n - i] := -1; end for;
      // Expand the vector of mult. of each generator & recompute.
      E := [InsertBlock(ZeroMatrix(ZZ, 1, n), e, 1, 1): e in E];
      for i in [1..k] do E[g][1][n - i + 1] := E[g][1][p]; end for;
      V := [e*Transpose(P^-1) : e in E]; v := ZeroMatrix(ZZ, 1, n);
      for i in [1..n] do v[1][i] := Min([vj[1][i] : vj in Prune(V)]); end for;
    end if;
  end for;
  // ------------ Add NEW satellite points ------------------
  e := v*Transpose(P); inCluster := [i : i in [1..n] | e[1][i] ne 0];
  points2test := #inCluster - 1; i := 2;
  while points2test gt 0 do
    // Values for the generators at point p.
    p := inCluster[i]; Vp := [vi[1][p] - v[1][p] : vi in Prune(V)];
    // Points p is proximate to && Points proximate to p.
    p_prox := [i : i in [1..#inCluster] | P[p][inCluster[i]] eq -1];
    prox_p := [i : i in [1..#inCluster] | P[inCluster[i]][p] eq -1];
    Q := [ q : q in p_prox | &+Eltseq(Submatrix(P, prox_p, [q])) eq 0];
    for q in Q do
      // Values for the generators at point q.
      Vq := [vi[1][q] - v[1][q] : vi in Prune(V)];
      if &*[Vp[i] + Vq[i] : i in [1..#Vp]] ne 0 then
        // Expand the proximity matrix.
        P := InsertBlock(ScalarMatrix(n + 1, 1), P, 1, 1); n := n + 1;
        P[n][p] := -1; P[n][q] := -1;
        // Expand the vector of mult. of each generator.
        E := [InsertBlock(ZeroMatrix(ZZ, 1, n), e, 1, 1): e in E];
        V := [e*Transpose(P^-1) : e in E]; v := ZeroMatrix(ZZ, 1, n);
        for i in [1..n] do v[1][i] := Min([vj[1][i] : vj in Prune(V)]); end for;
        e := v*Transpose(P); inCluster := [i : i in [1..n] | e[1][i] ne 0];
        points2test := points2test + 1;
      end if;
    end for;
    points2test := points2test - 1; i := i + 1;
  end while;
  return <P, v + V[#V]>;
end intrinsic;
