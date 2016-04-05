
Euclides := function(m, n)
  hs := []; ns := [];
  while n ne 0 do
    Append(~hs, m div n); Append(~ns, n);
    r := m mod n; m := n; n := r;
  end while; return <hs, ns>;
end function;

intrinsic CharExponents(s::RngSerPuisElt) -> []
{ Returns the characteristic exponents of a Puiseux series }
  C, m, n := ElementToSequence(s);
  // Exponents appearing in the series s
  E := [ m + i - 1 : i in [1 .. #C] | C[i] ne 0 ];
  charExps := [<0, n>]; ni := n;
  while ni ne 1 do
    // m_i = min{ j | a_j != 0 and j \not\in (n_{i-1}) }
    mi := [ e : e in E | e mod ni ne 0 ][1];
    // n_i = gcd(n, m_1, ..., m_k)
    ni := Gcd(ni, mi);
    Append(~charExps, <mi, ni>);
  end while; return charExps;
end intrinsic;

intrinsic CharExponents(f::RngMPolElt) -> []
{ Returns the characteristic exponents of an irreducible bivariate polynomials }
require Rank(Parent(f)) eq 2: "Argument must be a bivariate polynomial";
  S := NewtonPuiseuxAlgorithmReduced(f);
  if #S ne 1 then error "Argument must be an irreducible series"; end if;
  return CharExponents(S[1]);
end intrinsic;

intrinsic CharExponents(G::[ RngIntElt ]) -> [ ]
{ Computes the characteristic exponents from the generators of the semigroup }
  M := [ G[1] ]; N := [ G[1] ];
  for i in [2..#G] do
    M cat:= [ &+[j ne i select -(N[j - 1] - N[j]) div N[i - 1] * M[j]
      else G[j] : j in [2..i]] ]; N cat:= [ Gcd(M) ];
  end for;
  return [<0, N[1]>] cat [<M[i], N[i]> : i in [2..#M]];
end intrinsic;

TailExponentSeries := function(s)
  C, m, n := ElementToSequence(s);
  // Exponents appearing in the series s
  E := [ m + i - 1 : i in [1 .. #C] | C[i] ne 0 ];
  charExps := CharExponents(s); g := #charExps;
  if s eq 0 then return <0, 1>; end if;
  return [<e, 1> : e in E | e ge charExps[g][1]][1];
end function;

TailExponentMatrix := function(P, v)
  E := CharExponents(SemiGroup(P, v)); n := Ncols(P);
  // If last point is satellite there is no tail exponent
  if &+Eltseq(P[n]) eq -1 then return -1; end if;
  isSat := &+[Transpose(P)[j]: j in [1..n]];
  p := ([i : i in Reverse([1..n]) | isSat[i] eq -1] cat [0])[1] + 1;
  return E[#E][1] + (n - p);
end function;

intrinsic SemiGroup(s::RngSerPuisElt) -> []
{ Computes a minimal set of generators for the semigroup of the
  Puiseux series of an irreducible plane curve }
  M := CharExponents(s);
  // (G)amma starts with <n, m, ...>
  G := [i gt 2 select ( (Self(i-1) - M[i-1][1]) * M[i-2][2] div M[i-1][2] ) +
    M[i][1] + ( (M[i-2][2] - M[i-1][2]) div M[i-1][2] ) * M[i-1][1]
        else [M[1][2], M[2][1]][i] : i in [1 .. #M]];
  return G;
end intrinsic;

intrinsic SemiGroup(f::RngMPolElt) -> []
{ Computes a minimal set of generators for the semigroup of
  and irreducible plane curve }
require Rank(Parent(f)) eq 2: "Argument must be a bivariate polynomial";
  S := NewtonPuiseuxAlgorithmReduced(f);
  if #S ne 1 then error "Argument must be an irreducible series"; end if;
  return SemiGroup(S[1]);
end intrinsic;

intrinsic SemiGroup(P::Mtrx, v::Mtrx) -> []
{ Returns the minimal set of generators for the semigroup of and irreducible
  plane curve from its weighted cluster of singular points }
require Ncols(P) eq Nrows(P) and Ncols(v) eq Ncols(P) and Nrows(v) eq 1:
  "Arguments do not have the required dimensions";
require v[1][1] le v[1][Ncols(v)]: "Argument v is not a vector of values";
require #[i: i in [1..Ncols(P)] | (v*P)[1][i] gt 0] eq 1 or Gcd(Eltseq(v)) eq 1:
  "Weighted cluster not irreducible";

  G := [v[1][1]]; n := Ncols(P); isSat := &+[Transpose(P)[i] : i in [1..n]];
  G cat:= [v[1][i] : i in [1..n - 1] | isSat[i] ne -1 and isSat[i + 1] eq -1];
  return G;
end intrinsic;

SemiGroupCoord := function(v, G)
  return [];
end function;
