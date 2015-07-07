
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
----------- These are methods that M2 should have, but it doesn't. -------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

polyDivision = method(TypicalValue => Sequence)
polyDivision (RingElement, RingElement) := (u, v) -> (
  eps := 2.0^(-(precision ring u)/2);
  x := first generators ring u;
  m := (degree u)#0; n := (degree v)#0;
  q := 0;
  for kk from 0 to m - n do (
    k := m - n - kk;
    qk := u_(x^(n+k))//v_(x^n);
    q = q + qk*x^k;
    u = u - qk*x^k*v;
  ); return (clean(eps, q), clean(eps, u));
)

quot = method(TypicalValue => RingElement)
quot (RingElement, RingElement) := (u, v) -> first polyDivision(u, v);

rem = method(TypicalValue => RingElement)
rem (RingElement, RingElement) := (u, v) -> last polyDivision(u, v);

GCD = method(TypicalValue => RingElement);
GCD (RingElement, RingElement) := (a, b) -> (
  while b != 0 do (
    r := rem(a,b); a = b; b = r;
  ); return a / leadCoefficient(a);
)

squareFreePart = method(TypicalValue => RingElement);
squareFreePart (RingElement) := (f) -> (
  y := last generators ring f;
  return f//gcd(f, diff(y, f));
)

squareFreeFactorization = method(TypicalValue => List);
squareFreeFactorization (RingElement) := (f) -> (
  -- Yun's algorithm
  x := first generators ring f;
  y := last generators ring f;
  d := diff(y, f);
  squareFree := {};
  while degree f != {0} do (
    a := gcd(f, d);
    f = f//a;
    d = d//a - diff(y, f);
    squareFree = append(squareFree, a);
  );
  squareFree = drop(squareFree, 1);
  return apply(select(pack(2, mingle(squareFree, 1..#squareFree)),
    fact -> not isConstant fact#0), toSequence);
)

ZZ * InfiniteNumber := (n, inf) -> if n == 0 then 0 else
  if n > 0 then infinity else -infinity;

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

