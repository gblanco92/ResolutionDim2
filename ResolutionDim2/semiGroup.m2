
charExponents = method(TypicalValue => List);
charExponents (PuiseuxSeries) := (s) -> (
  exps := reverse exponents s.p;
  charExps := {(0, s.n)};
  ni := s.n;
  while ni != 1 do (
    -- m_i = min{ j | a_j != 0 and j \not\in (n_{i-1}) }
    idx := position(exps, e -> e#0 % ni != 0);
    mi := exps#idx#0;
    -- n_i = gcd(n, m_1, ..., m_k)
    ni = gcd(ni, mi);
    charExps = append(charExps, (mi, ni));
  ); return charExps;
)

charExponents (RingElement) := (f) -> (
  if not isPolynomialRing ring f then error "not a polynomial";
  if numgens ring f != 2 then error "not a bivariate polynomial";
  S := puiseuxExpansion(f);
  if #S > 1 then error "not an irreducible curve";
  (s, m) := first S;
  return charExponents s;
)

tailExponents = method(TypicalValue => List);
tailExponents (PuiseuxSeries) := (s) -> (
  exps := reverse exponents s.p;
  (mk, nk) := last charExponents(s);
  -- If smooth branch take all, otherwise look for the last char. exp.
  if mk == 0 then idx := 0 else idx = position(exps, e -> e#0 == mk);
  return take(apply(exps, e -> (e#0, 1)), {idx, #exps});
)

semiGroup = method(TypicalValue => List);
semiGroup (PuiseuxSeries) := (s) -> (
  C := charExponents s;
  -- (G)amma starts with <n, m, ...>
  G := {C#0#1, C#1#0};
  for i from 2 to #C-1 do (
    G = append(G, ( (G#(i-1) - C#(i-1)#0)*C#(i-2)#1/C#(i-1)#1 ) +
          C#i#0 + ( (C#(i-2)#1 - C#(i-1)#1)/C#(i-1)#1 )*C#(i-1)#0 );
  );
  return G;
)

semiGroup (RingElement) := (f) -> (
  if not isPolynomialRing ring f then error "not a polynomial";
  if numgens ring f != 2 then error "not a bivariate polynomial";
  S := puiseuxExpansion(f);
  if #S > 1 then error "not an irreducible curve";
  (s, m) := first S;
  return semiGroup s;
)

