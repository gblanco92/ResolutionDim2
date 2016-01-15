needs "./examples.m2"

S = puiseuxExpansion(book);
assert(#S == 4);
assert(listForm first S#0 == {({1, 0}, -1)});
assert(listForm first S#1 == {({1, 0}, 1)});
assert(listForm first S#2 == {({2, 0}, 1)});
assert(listForm first S#3 == {({3, 0}, -1)});

S = puiseuxExpansion(aa, Terms => 10);
assert(#S == 1);
assert(#listForm first S#0 == 2);
assert(last S#0 == 1);
assert(charExponents first S#0 == {(0, 4), (6, 2), (17, 1)});

S = puiseuxExpansion(bb, Terms => 5);
assert(#S == 1);
assert(#listForm first S#0 == 5);
assert(last S#0 == 1);
assert(charExponents first S#0 == {(0, 10), (4, 2), (5, 1)});

S = puiseuxExpansion(bbb^2);
assert(#S == 1);
assert(#listForm first S#0 == 2);
assert(last S#0 == 2);
assert(charExponents first S#0 == {(0, 4), (10, 2), (11, 1)});

S = puiseuxExpansion(cc^5);
assert(#S == 1);
assert(#listForm first S#0 == 2);
assert(last S#0 == 5);
assert(charExponents first S#0 == {(0, 4), (6, 2), (15, 1)});

S = puiseuxExpansion(aa*bb*cc*bbb);
assert(#S == 4);
assert(#listForm first S#0 == 2);
assert(#listForm first S#1 == 2);
assert(#listForm first S#2 == 2);
assert(#listForm first S#3 == 2);
assert(last S#0 == 1);
assert(last S#1 == 1);
assert(last S#2 == 1);
assert(last S#3 == 1);
assert(charExponents first S#0 == {(0, 10), (4, 2), (5, 1)});
assert(charExponents first S#1 == {(0, 4), (6, 2), (15, 1)});
assert(charExponents first S#2 == {(0, 4), (6, 2), (17, 1)});
assert(charExponents first S#3 == {(0, 4), (10, 2), (11, 1)});

S = puiseuxExpansion(aa^2*bb^3*cc^4);
assert(#S == 3);
assert(#listForm first S#0 == 2);
assert(#listForm first S#1 == 2);
assert(#listForm first S#2 == 2);
assert(last S#0 == 3);
assert(last S#1== 4);
assert(last S#2 == 2);
assert(charExponents first S#0 == {(0, 10), (4, 2), (5, 1)});
assert(charExponents first S#1 == {(0, 4), (6, 2), (15, 1)});
assert(charExponents first S#2 == {(0, 4), (6, 2), (17, 1)});

S = puiseuxExpansion(ee^2);
assert(#S == 2);
assert(#listForm first S#0 == 10);
assert(#listForm first S#1 == 10);
assert(last S#0 == 2);
assert(last S#1 == 2);
assert(charExponents first S#0 == {(0, 10), (4, 2), (5, 1)});
assert(charExponents first S#1 == {(0, 10), (4, 2), (5, 1)});

S = puiseuxExpansion(ff, Terms => 4);
assert(#S == 2);
assert(#listForm first S#0 == 4);
assert(#listForm first S#1 == 4);
assert(last S#0 == 1);
assert(last S#1 == 1);
assert(charExponents first S#0 == {(0, 4), (10, 2), (11, 1)});
assert(charExponents first S#1 == {(0, 4), (10, 2), (11, 1)});

S = puiseuxExpansion(gg);
assert(#S == 2);
assert(#listForm first S#0 == 15);
assert(#listForm first S#1 == 15);
assert(last S#0 == 1);
assert(last S#1 == 1);
assert(charExponents first S#0 == {(0, 10), (4, 2), (5, 1)});
assert(charExponents first S#1 == {(0, 10), (4, 2), (5, 1)});

S = puiseuxExpansion(mm);
assert(#S == 1);
assert(#listForm first S#0 == 4);
assert(last S#0 == 1);
assert(charExponents first S#0 == {(0, 30), (12, 6), (15, 3), (19, 1)});

S = puiseuxExpansion(x^3*ll);
assert(#S == 2);
assert(#listForm first S#0 == 1);
assert(#listForm first S#1 == 3);
assert(last S#0 == 3);
assert(last S#1 == 1);
assert(charExponents first S#1 == {(0, 12), (30, 6), (33, 3), (37, 1)});

S = puiseuxExpansion(ee*gg);
assert(#S == 4);
assert(#listForm first S#0 == 15);
assert(#listForm first S#1 == 15);
assert(#listForm first S#2 == 10);
assert(#listForm first S#3 == 10);
assert(last S#0 == 1);
assert(last S#1 == 1);
assert(last S#2 == 1);
assert(last S#3 == 1);
assert(charExponents first S#0 == {(0, 10), (4, 2), (5, 1)});
assert(charExponents first S#1 == {(0, 10), (4, 2), (5, 1)});
assert(charExponents first S#2 == {(0, 10), (4, 2), (5, 1)});
assert(charExponents first S#3 == {(0, 10), (4, 2), (5, 1)});

S = puiseuxExpansion(ff^2*jj^3);
assert(#S == 4);
assert(#listForm first S#0 == 3);
assert(#listForm first S#1 == 3);
assert(#listForm first S#2 == 3);
assert(#listForm first S#3 == 3);
assert(last S#0 == 2);
assert(last S#1 == 2);
assert(last S#2 == 3);
assert(last S#3 == 3);
assert(charExponents first S#0 == {(0, 4), (10, 2), (11, 1)});
assert(charExponents first S#1 == {(0, 4), (10, 2), (11, 1)});
assert(charExponents first S#2 == {(0, 4), (10, 2), (11, 1)});
assert(charExponents first S#3 == {(0, 4), (10, 2), (11, 1)});

S = puiseuxExpansion(f^4, Terms => 10);
assert(#S == 3);
assert(#listForm first S#0 == 1);
assert(#listForm first S#1 == 1);
assert(#listForm first S#2 == 1);
assert(last S#0 == 4);
assert(last S#1 == 4);
assert(last S#2 == 4);
assert(charExponents first S#0 == {(0, 1)});
assert(charExponents first S#1 == {(0, 1)});
assert(charExponents first S#2 == {(0, 1)});

S = puiseuxExpansion(h^2);
assert(#S == 5);
assert(#listForm first S#0 == 1);
assert(#listForm first S#1 == 1);
assert(#listForm first S#2 == 1);
assert(#listForm first S#3 == 1);
assert(#listForm first S#4 == 1);
assert(last S#0 == 2);
assert(last S#1 == 2);
assert(last S#2 == 2);
assert(last S#3 == 2);
assert(last S#4 == 2);
assert(charExponents first S#0 == {(0, 1)});
assert(charExponents first S#1 == {(0, 1)});
assert(charExponents first S#2 == {(0, 1)});
assert(charExponents first S#3 == {(0, 1)});
assert(charExponents first S#4 == {(0, 1)});

S = puiseuxExpansion(foo);
assert(#S == 1);
assert(#listForm first S#0 == 1);
assert(last S#0 == 1);
assert(charExponents first S#0 == {(0, 1)});

S = puiseuxExpansion(bar);
assert(#S == 2);
assert(#listForm first S#0 == 0);
assert(#listForm first S#1 == 1);
assert(last S#0 == 1);
assert(last S#1 == 1);
assert(charExponents first S#0 == {(0, 1)});
assert(charExponents first S#1 == {(0, 1)});
