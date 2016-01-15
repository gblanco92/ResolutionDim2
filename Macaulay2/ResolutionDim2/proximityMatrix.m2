proximityMatrix = method(TypicalValue => Sequence,
                         Options => { ExtraPoint => false, Bits => 300 });
proximityMatrix (RingElement) := opts -> (f) -> (
  if not isPolynomialRing ring f then error "not a polynomial";
  if numgens ring f != 2 then error "not a bivariate polynomial";
  -- Get the Puiseux expansion of f.
  branches := puiseuxExpansion(f, Bits => opts.Bits);
  return proximityMatrix(branches,
    ExtraPoint => opts.ExtraPoint, Bits => opts.Bits);
)

proximityMatrix (List) := opts -> (branches) -> (
  -- Compute the proximity matrix and the contact matrix of each branch.
  contactMat := contactMatrix(branches);
  -- Proximity matrix of each branch.
  branchProx := apply(branches, 0..#branches - 1, (s, i) ->
    proximityMatrixBranch(s#0, max flatten entries contactMat^{i},
      ExtraPoint => opts.ExtraPoint));
  -- Compute the multiplicities of the infinitely near points of each branch.
  branchMult := apply(branches, 0..#branches - 1, (s, i) -> s#1 *
    multiplicityVectorBranch(s#0, max flatten entries contactMat^{i},
      ExtraPoint => opts.ExtraPoint));
  -- Get the proximity matrix of f and the position of each infinitely
  -- near point inside P.
  (P, p) := proximityMatrix(contactMat, branchProx);
  mult := {};
  -- Rearranges each point's multiplicity so its position is coherent with P.
  for i from 0 to #branches - 1 do (
    m := new MutableList from (numcols(P):0);
    for j from 0 to #p#i - 1 do m#(p#i#j) = branchMult#i_j;
    mult = mult | {vector toList m};
  ); return (P, mult);
)

proximityMatrix (Matrix, List) := opts -> (contactMat, branchProx) -> (
  ------------------------- Base case --------------------------------
  -- If there is only branch, return its prox. matrix.
  if #branchProx == 1 then
    return (branchProx#0, {toList(0..numcols(branchProx#0)-1)});
  ------------------- Compute the splits -----------------------------
  -- Substract one to all the contact numbers except the diagonal ones.
  contactMat = contactMat - matrix pack(#branchProx, (#branchProx)^2:1) +
    matrix mutableIdentity(ZZ, #branchProx);
  -- Idenitify each current branch with an ID from 0 to #brances.
  C := contactMat;
  remainingBranch := toList(0..numcols(C) - 1);
  -- Splits will contain lists of branches ID, where two branches will
  -- be in the same list iff they don't separate in the current node.
  splits := {};
  while #remainingBranch != 0 do (
    -- Get the contact number of the first remaining branch;
    branchContacts := first entries C;
    -- Get the positions of the branches with contact > 1 & contact = 1.
    sameBranchIndex := positions(branchContacts, c -> c != 0);
    otherBranchIndex := positions(branchContacts, c -> c == 0);
    -- Save the branches with contact > 1 together.
    splits = append(splits, remainingBranch_sameBranchIndex);
    -- Remove those branches since they've been splitted from the rest.
    remainingBranch = remainingBranch_otherBranchIndex;
    -- Compute the contact matrix of the remaining branches.
    C = submatrix(C, otherBranchIndex, otherBranchIndex);
  );
  ---------- Compute the prox. matrix of each subdiagram -------------
  -- Substract one to all the contact numbers and erase the
  -- first point of the proximity matricies of the current
  -- branches since we are moving down the Enriques diagram.
  newBranchProx := apply(branchProx, P -> submatrix'(P, {0}, {0}));
  -- Traverse each sub-diagram recursively.
  splitResult := apply(splits, split -> (
    proximityMatrix(submatrix(contactMat, split, split), newBranchProx_split)));
  -------------- Merge the prox. matrix of each split ----------------
  -- Create the matrix that will hold the proximity branch of this subdiagram.
  numPoints := sum apply(splitResult, (M, pos) -> numcols M) + 1;
  P := mutableIdentity(ZZ, numPoints);
  rowPoint := {}; k := 0;
  -- For each set of branches that splits in this node...
  for s from 0 to #splits - 1 do (
    -- Get the proximity matrix & the position of the points
    -- (relative to that prox. matrix) of the s-th subdiagram.
    (M, splitRowPoint) := splitResult#s;
    -- Copy the submatrix M inside P with the top left entry in (k+1, k+1)
    copySubmatrix(P, M, k + 1);
    -- Sum k+1 and add the new point ({0}) to the position of the
    -- points relative to the prox. matrix of the subdiagram.
    splitRowPoint = apply(splitRowPoint, pp -> {0} | apply(pp, p -> p + k + 1));
    rowPoint = rowPoint | splitRowPoint;
    -- Use the information in splitRowPoint to set the proximities of
    -- the current point into the new prox. matrix (P):
    -- For each branch in this subdiagram...
    for i from 0 to #(splits#s) - 1 do (
      Q := branchProx#(splits#s#i);
      -- For each element int the first column...
      for j from 1 to numcols(Q) - 1 do P_((splitRowPoint#i)#j, 0) = Q_(j, 0);
    ); k = k + numcols(M);
  );
  -- Make sure rowPoint is returned in the original order.
  splits = flatten splits;
  splits = toList apply(0..#splits-1, i -> position(splits, j -> j == i));
  return (matrix P, rowPoint_splits);
)

contactMatrix = method(TypicalValue => Matrix);
contactMatrix (List) := (branches) -> (
  -- Add a dummy term so compare exact branches is easier.
  x := first generators ring (first branches)#0;
  maxExp := max apply(branches, (s, m) -> ceiling((last exponents s)#0 + 1));
  branches = apply(branches, (s, m) -> (s + x^maxExp, m));
  branchesInfo := apply(branches, (s, m) -> puiseuxInfo s);
  contact := mutableIdentity(ZZ, #branches);
  -- For each pair of branches compute their contact number.
  for i from 0 to #branches - 1 do (
    for j from i + 1 to #branches - 1 do (
      contactNum := contactNumber(branchesInfo#i, branchesInfo#j);
      contact_(i,j) = contact_(j,i) = contactNum;
    );
  ); return matrix contact;
)

contactNumber = method(TypicalValue => ZZ);
contactNumber (List, List) := (branchInfoA, branchInfoB) -> (
  contactNumber := 0;
  -- For each characteristic exponent...
  for r from 0 to min(#branchInfoA, #branchInfoB) - 1 do (
    -- Get the contact number of this char. exponent and whether
    -- or not we should compare more points.
    (numExp, compNext) := contactNumberExp(branchInfoA#r, branchInfoB#r);
    contactNumber = contactNumber + numExp;
    if not compNext then break;
  ); return contactNumber;
)

contactNumberExp = method(TypicalValue => Sequence);
contactNumberExp (Sequence, Sequence) := (expInfoA, expInfoB) -> (
  contactNum := 0;
  -- Free points associated with the char. exponent.
  freeA := expInfoA#0;
  freeB := expInfoB#0;
  -- Satellite points associated with the char. exponent.
  satelliteA := new MutableList from expInfoA#1;
  satelliteB := new MutableList from expInfoB#1;
  -- Compare free points.
  for i from 0 to min(#freeA, #freeB) - 1 do (
    if freeA#i == freeB#i then contactNum = contactNum + 1
    else return (contactNum, false);
  );
  -- If the number of free points is not the same, no more points can be shared.
  if #freeA != #freeB then return (contactNum, false);
  -- Compare satellite points.
  satelliteA#-1 = satelliteA#-1 - 1;
  satelliteB#-1 = satelliteB#-1 - 1;
  for i from 1 to min(#satelliteA, #satelliteB) - 1 do (
    contactNum = contactNum + min(satelliteA#i, satelliteB#i);
    if satelliteA#i != satelliteB#i then return (contactNum, false);
  );
  -- If the number of stairs is not the same, no more points can be shared.
  if #satelliteA != #satelliteB then return (contactNum, false);
  -- Otherwise, all the points are shared.
  return (contactNum, true);
)

puiseuxInfo = method(TypicalValue => List);
puiseuxInfo (PuiseuxSeries) := (s) -> (
  pInfo := {}; allExps := charExponents(s);
  if tailExponents(s) != {} then allExps = allExps | { last tailExponents(s) };
  x := first generators ring s;
  for i from 1 to #allExps - 1 do (
    (mj, nj) := allExps#(i-1);
    mi := allExps#i#0;
    h0 := (mi - mj)//nj;
    free := apply(apply(0..h0, l -> (mj + l*nj)/s.n), e -> (e, s_(x^e)));
    satellite := first euclides(mi - mj, nj);
    pInfo = append(pInfo, (toList free, satellite));
  ); return pInfo;
)

puiseuxInfo (RingElement) := (x) -> {({(0,0)}, {0, infinity})}

multiplicityVectorBranch = method(TypicalValue => List,
                                  Options => { ExtraPoint => false });
multiplicityVectorBranch (PuiseuxSeries, ZZ) := opts -> (s, maxContact) -> (
  mult := (); charExps := charExponents(s);
  for i from 1 to #charExps - 1 do (
    (mj, nj) := charExps#(i-1);
    mi := charExps#i#0;
    (hs, ns) := euclides(mi - mj, nj);
    scan(hs, ns, (h,n) -> mult = mult | (h:n));
  );
  mult = mult | ((maxContact - #mult):1);
  if opts.ExtraPoint then mult = mult | (1:1);
  return vector toList mult;
)

multiplicityVectorBranch (RingElement, ZZ) := opts -> (x, maxContact) -> (
  if opts.ExtraPoint then maxContact = maxContact + 1;
  return vector toList (maxContact:1);
)

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
    charExps = append(charExps, (exps#idx#0, ni));
  ); return charExps;
)

tailExponents = method(TypicalValue => List);
tailExponents (PuiseuxSeries) := (s) -> (
  exps := reverse exponents s.p;
  (mk, nk) := last charExponents(s);
  -- If smooth branch take all, otherwise look for the last char. exp.
  if mk == 0 then idx := 0 else idx = position(exps, e -> e#0 == mk);
  return take(apply(exps, e -> (e#0, 1)), {idx, #exps});
)

euclides = method(TypicalValue => List);
euclides (ZZ, ZZ) := (m, n) -> (
  hs := {}; ns := {};
  while n != 0 do (
    hs = hs | {m // n};
    ns = ns | {n};
    r := m % n; m = n; n = r;
  ); return (hs, ns);
)

proximityMatrixBranch = method(TypicalValue => Matrix,
                               Options => {ExtraPoint => false});
proximityMatrixBranch (PuiseuxSeries, ZZ) := opts -> (branch, maxContact) -> (
  h := drop(apply(puiseuxInfo branch, charExps -> charExps#1), -1);
  numPoints := max(sum flatten h, maxContact);
  if opts.ExtraPoint then numPoints = numPoints + 1;
  -- Construct a prioximity matrix with free points only.
  prox := expandMatrix(matrix {{}}, numPoints);
  -- Fill in satellite points proximities.
  for i from 0 to #h - 1 do (
    -- Inverted branch case.
    if i == 0 and h#0#0 == 0 then start := 2 else start = 1;
    hi := new MutableList from h#i; hi#-1 = hi#-1 - 1;
    for j from start to #hi - 1 do (
      l := sum flatten h_{0..i-1} + sum (toList hi)_{0..j-1};
      for k from 1 to hi#j do prox_(l+k, l-1) = -1;
    );
  ); return matrix prox;
)

proximityMatrixBranch (RingElement, ZZ) := opts -> (branch, maxContact) -> (
  if opts.ExtraPoint then maxContact = maxContact + 1;
  return matrix expandMatrix(matrix {{}}, maxContact);
)

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------

copySubmatrix = method(TypicalValue => Nothing);
copySubmatrix (MutableMatrix, Matrix, ZZ) := (A, B, k) -> (
  for i from 0 to numcols(B) - 1 do
    for j from 0 to numcols(B) - 1 do
      A_(k + i, k + j) = B_(i, j);
)

expandMatrix = method(TypicalValue => MutableMatrix)
expandMatrix (Matrix, ZZ) := (A, k) -> (
  newA := mutableIdentity(ZZ, numcols A + k);
  for i from 1 to numcols(newA) - 1 do newA_(i, i-1) = -1;
  copySubmatrix(newA, matrix A, 0);
  return newA;
)

expandVector = method(TypicalValue => Vector)
expandVector (Vector, ZZ) := (v, k) -> return vector(entries v | toList(k:0));
