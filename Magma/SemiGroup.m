
forward TailExponents;

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
  end while;
  return charExps;
end intrinsic;

intrinsic CharExponents(f::RngMPolElt) -> []
{ Returns the characteristic exponents of an irreducible bivariate polynomials }
require Rank(Parent(f)) eq 2: "Argument must be a bivariate polynomial";
  S := NewtonPuiseuxAlgorithmReduced(f);
  if #S ne 1 then error "Argument must be an irreducible series"; end if;
  return CharExponents(S[1]);
end intrinsic;

TailExponents := function(s)
  C, m, n := ElementToSequence(s);
  // Exponents appearing in the series s
  E := [ m + i - 1 : i in [1 .. #C] | C[i] ne 0 ];
  charExps := CharExponents(s); g := #charExps;
  if s eq 0 then return [<0, 1>]; end if;
  return [<e, 1> : e in E | e ge charExps[g][1]];
end function;

intrinsic SemiGroup(s::RngSerPuisElt) -> []
{ Computes a minimal set of generators for the semigroup of the
  Puiseux series of an irreducible plane curve }
  M := CharExponents(s);
  // (G)amma starts with <n, m, ...>
  G := [i gt 2 select ( (Self(i-1) - M[i-1][1]) * M[i-2][2] / M[i-1][2] ) +
    M[i][1] + ( (M[i-2][2] - M[i-1][2]) / M[i-1][2] ) * M[i-1][1]
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
