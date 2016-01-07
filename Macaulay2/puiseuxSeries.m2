--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
----------------------------- CONSTRUCTORS -------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

PuiseuxSerie = new Type of HashTable;

puiseuxSerie = method(TypicalValue => PuiseuxSerie);
puiseuxSerie (RingElement, ZZ) := (f, m) -> (
  if not isPolynomialRing ring f then error "not a polynomial";
  if (degree f)#0 <= 0 then return new PuiseuxSerie from hashTable {p=>f, n=>1};
  return new PuiseuxSerie from hashTable {p=>f, n=>m};
)

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
------------------------- OVERLOADED UNARY METHODS -----------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

ring (PuiseuxSerie) := (f) -> ring f.p;

terms (PuiseuxSerie) := (f) ->
  reverse apply(terms f.p, term -> puiseuxSerie(term, f.n));

PuiseuxSerie _ PuiseuxSerie := (f, m) -> (
  if size m != 1 then error "expected a monomial";
  R := coefficientRing ring f;
  e := first exponents m;
  idx := position(exponents f, i -> i == e);
  if idx === null then return 0_R else return (listForm f)#idx#1;
)

listForm (PuiseuxSerie) := (f) ->
  return reverse apply(listForm f.p, (e, c) -> (e/f.n, c));

exponents (PuiseuxSerie) := (f) -> if f.p == 0 then return {{0}} else
  reverse apply(exponents f.p, e -> {(first e)/f.n});

coefficients (PuiseuxSerie) := (f) -> coefficients f.p;

size (PuiseuxSerie) := (f) -> size f.p;

lift (PuiseuxSerie, PolynomialRing) := (f, R) -> (
  if f.n != 1 then error "cannot lift to polynomial ring";
  return lift(f.p, R);
)

- PuiseuxSerie := (f) -> puiseuxSerie(-f.p, f.n);

--rootUnity = method(TypicalValue => RingElement, Options => {Bits => 300})
--rootUnity (QQ) := opts -> (kn) -> (
--  bits := 2*opts.Bits;
--  ppi := numeric_bits pi;
--  root := exp(2*ppi*ii*kn);
--  if (abs(realPart root) < 2.0^(-bits/2)) then root = (imaginaryPart root)*ii;
--  if (abs(imaginaryPart root) < 2.0^(-bits/2)) then root = realPart root;
--  return root;
--)
--
--conjugate (PuiseuxSerie) := (f) -> (
--  x := first generators ring f;
--  n := f.n;
--  bits := precision coefficientRing ring f;
--  return apply(toList(1..n), k -> sum apply(listForm f, (e, c) ->
--    c*rootUnity((e#0)*k, Bits => bits//2)*x^(e#0)));
--)
--
--clean (PuiseuxSerie) := (f) -> (
--  x := first generators ring f; y := last generators ring f;
--  bits := precision coefficientRing ring f;
--  s := puiseuxSerie(0*x, 1);
--  s = s + sum apply(listForm f, (e, c) -> (
--    if (abs(realPart c) < 2.0^(-bits)) then c = (imaginaryPart c)*ii;
--    if (abs(imaginaryPart c) < 2.0^(-bits)) then c = realPart c;
--    return c*x^(e#0)*y^(e#1);
--  )); return s;
--)
--
--toPolynomial = method(TypicalValue => RingElement)
--toPolynomial (PuiseuxSerie) := (s) -> (
--  R := ring s;
--  y := puiseuxSerie(last generators R, 1);
--  factors := apply(conjugate s, si -> y - si);
--  p := clean(product apply(conjugate s, si -> y - si));
--  x := first generators R; y = last generators R;
--  return sum apply(listForm p.p, (e, a) -> a*x^(e#0//p.n)*y^(e#1//p.n));
--)

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
------------------------- OVERLOADED BINARY METHODS ----------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

PuiseuxSerie + PuiseuxSerie := (f, g) -> (
  n := lcm(f.n, g.n);
  fSubs := apply(generators ring f, gen -> gen => gen ^ (n // f.n));
  gSubs := apply(generators ring g, gen -> gen => gen ^ (n // g.n));
  ff := sub(f.p, fSubs);
  gg := sub(g.p, gSubs);
  return puiseuxSerie(ff + gg, n);
)

PuiseuxSerie - PuiseuxSerie := (f, g) -> f + (-g);

PuiseuxSerie * PuiseuxSerie := (f, g) -> (
  n := lcm(f.n, g.n);
  fSubs := apply(generators ring f, gen -> gen => gen ^ (n // f.n));
  gSubs := apply(generators ring g, gen -> gen => gen ^ (n // g.n));
  ff := sub(f.p, fSubs);
  gg := sub(g.p, gSubs);
  return puiseuxSerie(ff * gg, n);
)

Number + PuiseuxSerie := (n, f) -> puiseuxSerie(n + f.p, f.n);

PuiseuxSerie + Number := (f, n) -> n + f;

Number * PuiseuxSerie := (n, f) -> puiseuxSerie(n*f.p, f.n);

PuiseuxSerie * Number := (f, n) -> n * f;

RingElement + PuiseuxSerie := (p, f) -> puiseuxSerie(p, 1) + f;

PuiseuxSerie + RingElement := (f, p) -> f + puiseuxSerie(p, 1);

RingElement * PuiseuxSerie := (p, f) -> puiseuxSerie(p, 1) * f;

PuiseuxSerie * RingElement := (f, p) -> f * puiseuxSerie(p, 1);

PuiseuxSerie ^ ZZ := (f, k) -> puiseuxSerie(f.p ^ k, f.n);

RingElement ^ QQ := (f, q) -> (
  if size f > 1 then error "no method for binary operator ^ applied
    to polynomials with more than one monomial";
  return puiseuxSerie(f^(numerator q), denominator q);
)

PuiseuxSerie ^ QQ := (f, q) -> (
  if size f > 1 then error "no method for binary operator ^ applied
    to polynomials with more than one monomial";
  return puiseuxSerie(f.p^(numerator q), f.n*(denominator q));
)

substitute (PuiseuxSerie, Option) := (f, opt) -> (
  if not isPolynomialRing ring opt#1 or size opt#1 != 1 then
    error "substitution error";
  if class opt#1 === PuiseuxSerie then
    return puiseuxSerie(sub(f.p, opt#0 => (opt#1).p), f.n * (opt#1).n)
  else
    return puiseuxSerie(sub(f.p, opt), f.n);
)

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------- FORMATTING OUTPUT ----------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

expression PuiseuxSerie := f -> (
  gens := generators ring f;
  if length listForm f == 0 then return expression 0
  else return sum apply(listForm f,
      (exps, coef) -> coef * product apply(exps, gens, (e, g) ->
        (expression g)^e)) + O"("(expression first gens)^(1/f.n)")";

)

net PuiseuxSerie := f -> net expression f;

toString PuiseuxSerie := f -> toString expression f;

tex PuiseuxSerie := f -> tex expression f;

html PuiseuxSerie := f -> html expression f;

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
