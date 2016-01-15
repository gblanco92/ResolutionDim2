basePoints = method(TypicalValue => Sequence)
basePoints (Ideal) := (I) -> (
  gen := first entries generators I;
  if not all(gen, f -> numgens ring f == 2) then
    error "not bivariate polynomials";
  -- Remove the fixed part.
  if #gen > 1 then gen = gen//gcd(gen)
  else (
    (P, e) := proximityMatrix(gen#0);
    return (P, sum e, P^-1*(sum e));
  );
  -- Compute the Puiseux expansion & the prox. matrix for the product.
  branches := puiseuxExpansion(gen);
  (P, e) = proximityMatrix(apply(branches, (s, l) -> (s, 1)),
              ExtraPoint => true);
 -- Multiplicities of each generator.
  ee := new MutableList from apply(1..#gen, i -> vector toList(numcols P:0));
  for i from 0 to #branches - 1 do (
    (s, l) := branches#i;
    scan(l, (j, n) -> ee#(j-1) = ee#(j-1) + n * e#i);
  ); ee = toList ee;
  -- Values for each generator.
  vv := P^-1*ee;
  -- Values for the points in the cluster.
  v := vector apply(entries matrix vv, min);
  -- Multiplicities for the cluster base points.
  m := P*v;
  -- Remove points not in the cluster of base points.
  inCluster := positions(entries m, x -> x != 0);
  P = P_inCluster^inCluster;
  vv = apply(vv, vg -> vector apply(inCluster, i -> vg_i));
  ee = P*vv;
  v = vector apply(entries matrix vv, min);
  -- Add NEW free points.
  -- For each last free point on a branch...
  apply(positions(sum entries P, i -> i == 1), p -> (
    -- Values for each generator in the point p.
    vvp := (entries matrix vv)#p;
    uniqueGen := #select(vvp, x -> x == v_p) == 1;
    -- Index of the generator achieving the minimum.
    g := minPosition(vvp);
    -- Multiplicity of that generator in the point p.
    mgp := (entries matrix ee)#p#g;
    -- If there exist a unique generator achieving the minimum
    -- value vp and the excess is positive add new points...
    if uniqueGen and mgp != 0 then (
      -- Minimum of the values for all the generators but g.
      wp := min(delete(v_p, vvp));
      -- Number of new free points.
      k := ceiling((wp - v_p)/mgp);
      -- Expand the prox. matrix.
      n := numcols P;
      P = expandMatrix(P, k);
      P_(n, n - 1) = 0; P_(n, p) = -1;
      P = matrix P;
      -- Expand the vector of mult. of each gen. with zeros.
      ee = new MutableList from apply(ee, e -> expandVector(e, k));
      -- For the generator g fill the new points with mult. mp.
      eeg := new MutableList from entries ee#g;
      for j from 0 to k - 1 do eeg#(n + j) = mgp;
      ee#g = vector toList eeg;
      ee = toList ee;
      -- Recompute all the other vectors for the next iteration.
      vv = P^-1*ee;
      v = vector apply(entries matrix vv, min);
    );
  ));
  -- Add NEW satellite points.
  points2test := numcols P - 1; p := 1;
  while points2test != 0 do (
    -- Values for the generators at point p.
    vvp := apply((entries matrix vv)#p, x -> x - v_p);
    -- Points p is prox to. && Points prox. to p
    Pprox := positions(flatten entries P^{p}, x -> x == -1);
    proxP := positions(flatten entries P_{p}, x -> x == -1);
    apply(select(Pprox, q -> sum flatten entries P^proxP_{q} == 0), q -> (
      vvq := apply((entries matrix vv)#q, x -> x - v_q);
      if product(vvq + vvp) != 0 then (
        -- Expand proximity matrix with a new point.
        n := numcols P;
        P = expandMatrix(P, 1);
        P_(n, n - 1) = 0; P_(n, p) = -1; P_(n, q) = -1;
        P = matrix P;
        -- Expand de vector of multiplicities of each generator.
        ee = apply(ee, e -> expandVector(e, 1));
        -- Recompute all the other vector for the next iteration.
        vv = P^-1*ee;
        v = vector apply(entries matrix vv, min);
        points2test = points2test + 1;
      );)
    ); points2test = points2test - 1; p = p + 1;
  ); return (P, P*v, v);
)
