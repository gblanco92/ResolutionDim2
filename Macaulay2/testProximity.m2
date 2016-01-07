needs "./examples.m2"
needs "./proximityMatrix.m2"

(P, e) = proximityMatrix(x*y);
assert(numRows P == 1);
assert(numColumns P == 1);
assert(#e == 2);
assert(sum e == vector {2});
print ".";

(P, e) = proximityMatrix(x^2*y^3, ExtraPoint => true);
assert(numRows P == 3);
assert(numColumns P == 3);
assert(#e == 2);
assert(sum e == vector {5, 2, 3});
print ".";

(P, e) = proximityMatrix(x*y*(x-y)*(x^3 - y^2));
assert(numRows P == 3);
assert(numColumns P == 3);
assert(#e == 4);
assert(sum e == vector {5, 2, 1});
print ".";

(P, e) = proximityMatrix(x*y*(x-y)*(x^3 - y^2), ExtraPoint => true);
assert(numRows P == 7);
assert(numColumns P == 7);
assert(#e == 4);
assert(sum e == vector {5, 1, 2, 1, 1, 1, 1});
print ".";

(P, e) = proximityMatrix(x*aa*bb);
assert(numRows P == 15);
assert(numColumns P == 15);
assert(#e == 3);
assert(sum e == vector {9, 5, 3, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1});
print ".";

(P, e) = proximityMatrix(f1*f2*f3, ExtraPoint => true);
assert(numRows P == 10);
assert(numColumns P == 10);
assert(#e == 3);
assert(sum e == vector {8, 1, 1, 1, 3, 3, 2, 1, 1, 1});
print ".";

(P, e) = proximityMatrix(v1, ExtraPoint => true);
assert(numRows P == 5);
assert(numColumns P == 5);
assert(#e == 2);
assert(sum e == vector {4, 2, 2, 1, 1});
assert(P^-1*sum e == vector {4, 6, 12, 13, 13});
print ".";

(P, e) = proximityMatrix(v2);
assert(numRows P == 3);
assert(numColumns P == 3);
assert(#e == 2);
assert(sum e == vector {4, 2, 2});
assert(P^-1*sum e == vector {4, 6, 12});
print ".";

(P, e) = proximityMatrix(v3, ExtraPoint => true);
assert(numRows P == 7);
assert(numColumns P == 7);
assert(#e == 2);
assert(sum e == vector {4, 1, 1, 1, 1, 1, 1})
assert(P^-1*sum e == vector {4, 5, 10, 11, 5, 10, 11});
print ".";

(P, e) = proximityMatrix(v4);
assert(numRows P == 6);
assert(numColumns P == 6);
assert(#e == 1);
assert(sum e == vector {3, 3, 3, 2, 1, 1});
assert(P^-1*sum e == vector {3, 6, 9, 11, 21, 33});
print ".";

(P, e) = proximityMatrix(v5, ExtraPoint => true);
assert(numRows P == 7);
assert(numColumns P == 7);
assert(#e == 1);
assert(sum e == vector {3, 3, 3, 2, 1, 1, 1});
assert(P^-1*sum e == vector {3, 6, 9, 11, 21, 33, 34});
print ".";

(P, e) = proximityMatrix(v6);
assert(numRows P == 7);
assert(numColumns P == 7);
assert(#e == 2);
assert(sum e == vector {6, 2, 1, 1, 2, 1, 1})
assert(P^-1*sum e == vector {6, 8, 15, 24, 8, 15, 24});
print ".";

(P, e) = proximityMatrix(v7, ExtraPoint => true);
assert(numRows P == 10);
assert(numColumns P == 10);
assert(#e == 6);
assert(sum e == vector {6, 2, 1, 1, 2, 1, 1, 2, 1, 1});
assert(P^-1*sum e == vector {6, 8, 9, 9, 8, 9, 9, 8, 9, 9});
print ".";

(P, e) = proximityMatrix(v8);
assert(numRows P == 6);
assert(numColumns P == 6);
assert(#e == 1);
assert(sum e == vector {3, 3, 3, 1, 1, 1});
assert(P^-1*sum e == vector {3, 6, 9, 10, 20, 30});
print ".";

(P, e) = proximityMatrix(v9, ExtraPoint => true);
assert(numRows P == 5);
assert(numColumns P == 5);
assert(#e == 2);
assert(sum e == vector {3, 2, 1, 1, 1});
assert(P^-1*sum e == vector {3, 5, 6, 9, 10});
print ".";

(P, e) = proximityMatrix(v10);
assert(numRows P == 3);
assert(numColumns P == 3);
assert(#e == 3);
assert(sum e == vector {4, 3, 1});
assert(P^-1*sum e == vector {4, 7, 12});
print ".";

(P, e) = proximityMatrix(v11, ExtraPoint => true);
assert(numRows P == 8);
assert(numColumns P == 8);
assert(#e == 4);
assert(sum e == vector {5, 4, 2, 1, 1, 1, 1, 1});
assert(P^-1*sum e == vector {5, 9, 11, 12, 12, 15, 16, 10});
print ".";

(P, e) = proximityMatrix(v12);
assert(numRows P == 5);
assert(numColumns P == 5);
assert(#e == 2);
assert(sum e == vector {4, 3, 1, 1, 1});
assert(P^-1*sum e == vector {4, 7, 12, 8, 16});
print ".";

(P, e) = proximityMatrix(v13, ExtraPoint => true);
assert(numRows P == 8);
assert(numColumns P == 8);
assert(#e == 3);
assert(sum e == vector {5, 4, 2, 1, 1, 1, 1, 1});
assert(P^-1*sum e == vector {5, 9, 11, 12, 21, 22, 15, 16});
print ".";

(P, e) = proximityMatrix(v14);
assert(numRows P == 5);
assert(numColumns P == 5);
assert(#e == 2);
assert(sum e == vector {6, 4, 2, 2, 2});
assert(P^-1*sum e == vector {6, 10, 18, 30, 32});
print ".";

(P, e) = proximityMatrix(v15, ExtraPoint => true);
assert(numRows P == 5);
assert(numColumns P == 5);
assert(#e == 3);
assert(sum e == vector {3, 3, 1, 1, 1});
assert(P^-1*sum e == vector {3, 6, 7, 7, 7});
print ".";

(P, e) = proximityMatrix(v16);
assert(numRows P == 5);
assert(numColumns P == 5);
assert(#e == 5);
assert(sum e == vector {9, 7, 2, 2, 2});
assert(P^-1*sum e == vector {9, 16, 27, 45, 47});
print ".";

(P, e) = proximityMatrix(v17, ExtraPoint => true);
assert(numRows P == 12);
assert(numColumns P == 12);
assert(#e == 7);
assert(sum e == vector {28, 21, 7, 7, 7, 1, 1, 1, 1, 1, 1, 1});
assert(P^-1*sum e == vector {28, 49, 84, 140, 196, 197, 197, 197, 197, 197, 197, 197});
print ".";

(P, e) = proximityMatrix(v18);
assert(numRows P == 7);
assert(numColumns P == 7);
assert(#e == 1);
assert(sum e == vector {7, 7, 2, 2, 2, 1, 1});
assert(P^-1*sum e == vector {7, 14, 16, 32, 48, 63, 112});
print ".";

(P, e) = proximityMatrix(v19, ExtraPoint => true);
assert(numRows P == 9);
assert(numColumns P == 9);
assert(#e == 1);
assert(sum e == vector {10, 10, 3, 3, 3, 1, 1, 1, 1});
assert(P^-1*sum e == vector {10, 20, 23, 46, 69, 90, 160, 230, 231});
print ".";

(P, e) = proximityMatrix(v20);
assert(numRows P == 5);
assert(numColumns P == 5);
assert(#e == 3);
assert(sum e == vector {5, 1, 1, 1, 1});
assert(P^-1*sum e == vector {5, 6, 12, 6, 12});
print ".";

(P, e) = proximityMatrix(vA);
assert(numRows P == 12);
assert(numColumns P == 12);
assert(#e == 5);
assert(sum e == vector {32, 32, 15, 12, 4, 1, 3, 2, 1, 1, 2, 1});
assert(P^-1*sum e == vector {32, 64, 79, 155, 223, 288, 381, 538, 920, 694, 236, 316});
print ".";

(P, e) = proximityMatrix(vB);
assert(numRows P == 13);
assert(numColumns P == 13);
assert(#e == 4);
assert(sum e == vector {45, 33, 12, 11, 1, 10, 10, 5, 3, 1, 1, 1, 1});
assert(P^-1*sum e == vector {45, 78, 135, 224, 360, 312, 322, 327, 652, 975, 1628, 980, 1308});
print ".";

(P, e) = proximityMatrix(vC);
assert(numRows P == 15);
assert(numColumns P == 15);
assert(#e == 5);
assert(sum e == vector {50, 50, 32, 13, 3, 2, 1, 10, 9, 9, 6, 3, 2, 1, 1});
assert(P^-1*sum e == vector {50, 100, 132, 245, 348, 450, 799, 387, 528, 537, 543, 1083, 1628, 2712, 2172});
print ".";
