-- Power sums, Gorenstein algebras, and determinantal loci.
-- 1.3 Binary forms and Hankel Matrices
--- Example 1.34 (1)
R = QQ[X,Y];
F1 = X^3+4*X^2*Y+4*X*Y^2
Phi1 = inverseSystem F1

---- Length two decomposition
R = QQ[X,Y];
isSubset(ideal((2*X-Y)^2),sub(Phi1,R)) -- (2X-Y)^2= is an element of Ann(F1)
R11 = QQ[a,b][X,Y];
F1Unknown1 = (a*X+b*Y)*(X+2*Y)^2
diff(last flatten entries basis(3,R11),F1Unknown1)
factor sub(sub(F1Unknown1,{a=>1,b=>0}),R) -- found a GAD of F1 of length 2

---- Length three decomposision
R = QQ[X,Y];
isSubset(ideal(X^2*(4*X-3*Y)),sub(Phi1,R)) -- X^2(4X-3Y) is an element of Ann(F1)
R12 = QQ[a,b,c][X,Y];
F1Unknown2 = (a*X+b*Y)*Y^2+c*(3*X+4*Y)^3
DiffList = flatten entries diff(basis(2,R12),F1Unknown2);
DiffAnn = 4*DiffList#0-4*DiffList#1+DiffList#2
Sol = gens ker transpose last coefficients transpose sub(last coefficients DiffAnn,coefficientRing R12)
sub(1/27*sub(F1Unknown2,{a=>-36,b=>-64,c=>1}),R)

--- Example 1.34 (2), (3)
R = QQ[X,Y];
F2 = X^4+Y^4+(X+Y)^4
Phi2 = inverseSystem F2 -- Ann(F2)_2 is empty and generators of Ann(F2)_3 is not unique.
F3 = X^4+Y^4
Phi3 = inverseSystem F3 -- There is unique generator in Ann(F3)_2 and no further normalized GAD's of length 3.

--- Additional discussions
---- How to find apolar ideal? Example: F2
F2 = X^4+Y^4+(X+Y)^4;
mingens ker diff(basis(1,R), F2) -- no column vectors with degree 0 entries
mingens ker diff(basis(2,R), F2) -- no column vectors with degree 0 entries
mingens ker diff(basis(3,R), F2) -- two column vectors with degree 0 entries
(super basis(3,R)) * submatrix(mingens ker diff(basis(3,R), F2),{0,1}) -- Obtained the two generators of Ann(F)
mingens ker diff(basis(4,R), F2) -- four column vectors with degree 0 entries
(super basis(4,R)) * mingens ker diff(basis(4,R), F2) -- They form degree 4 polynomials that are apolar to F2. By adding F2, we obtain a basis of R_4.
mingens ker diff(basis(5,R), F2)
(super basis(5,R)) * mingens ker diff(basis(5,R), F2)

---- How to recover the homogeneous polynomial from its apolar ideal? Example: F2
ideal(super basis(4, inverseSystem(super basis(4,inverseSystem F2)))) == ideal F2
