import "ProximityMatrix.m": ProximityMatrixImpl, ProximityMatrixBranch,
                            MultiplicityVectorBranch, CoefficientsVectorBranch;
import "SemiGroup.m": TailExponentSeries;
import "IntegralClosure.m": IntegralClosureImpl, Unloading, CleanIdeal,
                            ClusterFactorization;

FiltrationImpl := function(s, f, M)
  // Compute the auto-intersection of the cluster of singular points K.
  e := ProximityMatrixImpl([<s, 1>])[2][1];
  KK := &+[ei*ei : ei in Eltseq(e)]; N := Ncols(e) + M - KK;
  // Get the proximity matrix with all the necessary points.
  s := NewtonPuiseuxAlgorithmExpandReduced(s, f: Terms := M - KK - 1)[1];
  P := ProximityMatrixBranch(s, N);
  e := MultiplicityVectorBranch(s, N);
  c := CoefficientsVectorBranch(s, N);
  // The first ideals of the filtrations are the maximal ideal.
  Q := PolynomialRing(CoefficientRing(Parent(s)), 2);
  R<x, y> := PolynomialRing(RationalField(), 2); n := e[1];
  H := [ideal<R | x, y> : i in [1..n]]; I := [1]; m_i := n;
  e_i := Matrix(1, 1, [1]); ZZ := IntegerRing();
  // Construct the i-th cluster.
  while m_i lt M do
    // Enlarge the previous cluster with one point on the curve.
    I cat:= [I[#I] + 1]; P_i := Submatrix(P, I, I);
    e_i := InsertBlock(ZeroMatrix(ZZ, 1, #I), e_i, 1, 1);
    e_i[1][#I] := 1; v_i := e_i*Transpose(P_i)^-1;
    // Unload K_i to get a strictly consistent cluster.
    v_i := Unloading(P_i, v_i); e_i := v_i*Transpose(P_i);
    I := [i : i in [1..Ncols(P_i)] | e_i[1][i] ne 0]; P_i := Submatrix(P_i, I, I);
    e_i := Submatrix(e_i, [1], I); v_i := Submatrix(v_i, [1], I); c_i := c[I];
    // Get the intersection [K, K_i] & and the complete ideal H_i.
    KK_i := &+[e_i[1][j]*e[j] : j in [1..#I]];
    HH_i := [ IntegralClosureImpl(K_j[1], K_j[2], K_j[3], Q) :
      K_j in ClusterFactorization(P_i, v_i, c_i) ];
    H_i := CleanIdeal([j gt 1 select Self(j - 1)*HH_i[j]
      else HH_i[1] : j in [1..#HH_i]][#HH_i]);
    // Fill the gaps in the filtration.
    H cat:= [ideal<R | Generators(H_i)> : i in [1..Min(KK_i, M) - m_i]];
    m_i := KK_i;
  end while; return H;
end function;

intrinsic Filtration(f::RngMPolElt, n::RngIntElt) -> []
{ Returns a filtration by complete ideals of an irreducible
  plane curve up to multiplicity n. }
require Rank(Parent(f)) eq 2: "First argument must be a bivariate polynomial";
require n ge 0: "Second argument must be a positive integer";

  S := NewtonPuiseuxAlgorithm(f: Polynomial := true);
  if #S gt 1 or S[1][2] gt 1 then error "the curve must be irreducible"; end if;
  s := S[1][1]; f := S[1][3]; e := ProximityMatrixImpl([<s, 1>])[2][1];
  KK := &+[ei*ei : ei in Eltseq(e)]; // Curve auto-intersection.

  return FiltrationImpl(s, f, n eq 0 select KK else n);
end intrinsic;

intrinsic JacobianFiltration(f::RngIntElt) -> []
{ Computes a filtration adapted to the Jacobian ideal
  of an irreducible plane curve. }
require Rank(Parent(f)) eq 2: "First argument must be a bivariate polynomial";

  S := NewtonPuiseuxAlgorithm(f: Polynomial := true);
  if #S gt 1 or S[1][2] gt 1 then error "the curve must be irreducible"; end if;
  s := S[1][1]; f := S[1][3];

end intrinsic;

