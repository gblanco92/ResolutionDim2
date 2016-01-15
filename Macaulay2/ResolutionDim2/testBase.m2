needs "./examples.m2"

jacobiIdeal = method(TypicalValue => Ideal);
jacobiIdeal (RingElement) := (f) -> (
  x := first generators ring f;
  y := last generators ring f;
  return ideal(diff(x, f), diff(y, f));
)

(P, e, v) = basePoints(jacobiIdeal v1);
assert(numRows P == 4);
assert(numColumns P == 4);
assert(#e == #v);
assert(e == vector {3, 2, 1, 1});

(P, e, v) = basePoints(v1 + jacobiIdeal v1);
assert(numRows P == 4);
assert(numColumns P == 4);
assert(#e == #v);
assert(e == vector {3, 2, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v2);
assert(numRows P == 4);
assert(numColumns P == 4);
assert(#e == #v);
assert(e == vector {3, 2, 1, 1});

(P, e, v) = basePoints(v2 + jacobiIdeal v2);
assert(numRows P == 4);
assert(numColumns P == 4);
assert(#e == #v);
assert(e == vector {3, 2, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v3);
assert(numRows P == 3);
assert(numColumns P == 3);
assert(#e == #v);
assert(e == vector {3, 1, 1});

(P, e, v) = basePoints(v3 + jacobiIdeal v3);
assert(numRows P == 3);
assert(numColumns P == 3);
assert(#e == #v);
assert(e == vector {3, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v4);
assert(numRows P == 8);
assert(numColumns P == 8);
assert(#e == #v);
assert(e == vector {2, 2, 2, 2, 1, 1, 1, 1});

(P, e, v) = basePoints(v4 + jacobiIdeal v4);
assert(numRows P == 8);
assert(numColumns P == 8);
assert(#e == #v);
assert(e == vector {2, 2, 2, 2, 1, 1, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v5);
assert(numRows P == 5);
assert(numColumns P == 5);
assert(#e == #v);
assert(e == vector {2, 2, 2, 2, 2});

(P, e, v) = basePoints(v5 + jacobiIdeal v5);
assert(numRows P == 5);
assert(numColumns P == 5);
assert(#e == #v);
assert(e == vector {2, 2, 2, 2, 2});

(P, e, v) = basePoints(jacobiIdeal v6);
assert(numRows P == 3);
assert(numColumns P == 3);
assert(#e == #v);
assert(e == vector {5, 2, 2});

(P, e, v) = basePoints(v6 + jacobiIdeal v6);
assert(numRows P == 3);
assert(numColumns P == 3);
assert(#e == #v);
assert(e == vector {5, 2, 2});

(P, e, v) = basePoints(jacobiIdeal v7);
assert(numRows P == 7);
assert(numColumns P == 7);
assert(#e == #v);
assert(e == vector {5, 1, 1, 1, 1, 1, 1});

(P, e, v) = basePoints(v7 + jacobiIdeal v7);
assert(numRows P == 7);
assert(numColumns P == 7);
assert(#e == #v);
assert(e == vector {5, 1, 1, 1, 1, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v8);
assert(numRows P == 6);
assert(numColumns P == 6);
assert(#e == #v);
assert(e == vector {2, 2, 2, 2, 1, 1});

(P, e, v) = basePoints(v8 + jacobiIdeal v8);
assert(numRows P == 6);
assert(numColumns P == 6);
assert(#e == #v);
assert(e == vector {2, 2, 2, 2, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v9);
assert(numRows P == 4);
assert(numColumns P == 4);
assert(#e == #v);
assert(e == vector {2, 1, 1, 1});

(P, e, v) = basePoints(v9 + jacobiIdeal v9);
assert(numRows P == 4);
assert(numColumns P == 4);
assert(#e == #v);
assert(e == vector {2, 1, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v10);
assert(numRows P == 5);
assert(numColumns P == 5);
assert(#e == #v);
assert(e == vector {3, 2, 1, 1, 1});

(P, e, v) = basePoints(v10 + jacobiIdeal v10);
assert(numRows P == 5);
assert(numColumns P == 5);
assert(#e == #v);
assert(e == vector {3, 2, 1, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v11);
assert(numRows P == 8);
assert(numColumns P == 8);
assert(#e == #v);
assert(e == vector {4, 3, 1, 1, 1, 1, 1, 1});

(P, e, v) = basePoints(v11 + jacobiIdeal v11);
assert(numRows P == 8);
assert(numColumns P == 8);
assert(#e == #v);
assert(e == vector {4, 3, 1, 1, 1, 1, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v12);
assert(numRows P == 6);
assert(numColumns P == 6);
assert(#e == #v);
assert(e == vector {3, 2, 1, 1, 1, 1});

(P, e, v) = basePoints(v12 + jacobiIdeal v12);
assert(numRows P == 6);
assert(numColumns P == 6);
assert(#e == #v);
assert(e == vector {3, 2, 1, 1, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v13);
assert(numRows P == 9);
assert(numColumns P == 9);
assert(#e == #v);
assert(e == vector {4, 3, 1, 1, 1, 1, 1, 1, 1});

(P, e, v) = basePoints(v13 + jacobiIdeal v13);
assert(numRows P == 9);
assert(numColumns P == 9);
assert(#e == #v);
assert(e == vector {4, 3, 1, 1, 1, 1, 1, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v14);
assert(numRows P == 8);
assert(numColumns P == 8);
assert(#e == #v);
assert(e == vector {5, 4, 1, 1, 1, 1, 1, 1});

(P, e, v) = basePoints(v14 + jacobiIdeal v14);
assert(numRows P == 8);
assert(numColumns P == 8);
assert(#e == #v);
assert(e == vector {5, 4, 1, 1, 1, 1, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v15);
assert(numRows P == 4);
assert(numColumns P == 4);
assert(#e == #v);
assert(e == vector {2, 2, 1, 1});

(P, e, v) = basePoints(v15 + jacobiIdeal v15);
assert(numRows P == 4);
assert(numColumns P == 4);
assert(#e == #v);
assert(e == vector {2, 2, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v16);
assert(numRows P == 12);
assert(numColumns P == 12);
assert(#e == #v);
assert(e == vector {8, 6, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1});

(P, e, v) = basePoints(v16 + jacobiIdeal v16);
assert(numRows P == 12);
assert(numColumns P == 12);
assert(#e == #v);
assert(e == vector {8, 6, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v17);
assert(numRows P == 23);
assert(numColumns P == 23);
assert(#e == #v);

(P, e, v) = basePoints(v17 + jacobiIdeal v17);
assert(numRows P == 23);
assert(numColumns P == 23);
assert(#e == #v);

(P, e, v) = basePoints(jacobiIdeal v18);
assert(numRows P == 4);
assert(numColumns P == 4);
assert(#e == #v);
assert(e == vector {6, 6, 3, 3});

(P, e, v) = basePoints(v18 + jacobiIdeal v18);
assert(numRows P == 4);
assert(numColumns P == 4);
assert(#e == #v);
assert(e == vector {6, 6, 3, 3});

(P, e, v) = basePoints(jacobiIdeal v19);
assert(numRows P == 8);
assert(numColumns P == 8);
assert(#e == #v);
assert(e == vector {9, 9, 4, 4, 1, 1, 1, 1});

(P, e, v) = basePoints(v19 + jacobiIdeal v19);
assert(numRows P == 8);
assert(numColumns P == 8);
assert(#e == #v);
assert(e == vector {9, 9, 4, 4, 1, 1, 1, 1});

(P, e, v) = basePoints(jacobiIdeal v20);
assert(numRows P == 3);
assert(numColumns P == 3);
assert(#e == #v);
assert(e == vector {4, 1, 1});

(P, e, v) = basePoints(v20 + jacobiIdeal v20);
assert(numRows P == 3);
assert(numColumns P == 3);
assert(#e == #v);
assert(e == vector {4, 1, 1});

(P, e, v) = basePoints(jacobiIdeal vA);
assert(numRows P == 54);
assert(numColumns P == 54);
assert(#e == #v);

(P, e, v) = basePoints(vA + jacobiIdeal vA);
assert(numRows P == 54);
assert(numColumns P == 54);
assert(#e == #v);

-- Too long to run!
-------------------
-- (P, e, v) = basePoints(jacobiIdeal vB);
-- assert(numRows P == 52);
-- assert(numColumns P == 52);
-- assert(#e == #v);
--
-- (P, e, v) = basePoints(vB + jacobiIdeal vB);
-- assert(numRows P == 52);
-- assert(numColumns P == 52);
-- assert(#e == #v);
--
-- (P, e, v) = basePoints(jacobiIdeal vC);
-- assert(numRows P == 103);
-- assert(numColumns P == 103);
-- assert(#e == #v);
--
-- (P, e, v) = basePoints(vC + jacobiIdeal vC);
-- assert(numRows P == 103);
-- assert(numColumns P == 103);
-- assert(#e == #v);
