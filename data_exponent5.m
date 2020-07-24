// In the file we keep the data for the computations of the Wordell-Weil group of genus 2 curves with positive rank


_<x> := PolynomialRing(Rationals());


// The case I

CI1 := HyperellipticCurve(3*x^5 + 35 * 2^1);
JI1 := Jacobian(CI1);
GensI1 := [JI1![x^2 - 1141/675*x + 1331/225, 6389/10125*x - 62374/3375]];

CI3 := HyperellipticCurve(3*x^5 + 35 * 2^3);
JI3 := Jacobian(CI3);
GensI3 := [JI3![x^2 + 28/25*x - 84/25, 537/125*x + 1414/125]];


CI5 := HyperellipticCurve(3*x^5 + 35 * 2^5);
JI5 := Jacobian(CI5);
GensI5 := [JI5![x^2 - 16*x + 121/4, 435/4*x - 2191/8],JI5![x - 3, -43],JI5![x + 2, -32]];



CI7 := HyperellipticCurve(3*x^5 + 35 * 2^7);
JI7 := Jacobian(CI7);
GensI7 := [JI7![x^2 + 20*x + 120, 76*x + 1360]];


CI9 := HyperellipticCurve(3*x^5 + 35 * 2^9);
JI9 := Jacobian(CI9);
GensI9 := [JI9![x^2 - 3088/27*x - 15712/27, 175580/81*x + 848864/81]];



// The case II

CII := HyperellipticCurve(15*x^5 + 35 * 2^5);
JII := Jacobian(CII);
GensII := [JII![x - 2, -40]];



// The case III

CIII1 := HyperellipticCurve(5*x^5 + 280);
JIII1 := Jacobian(CIII1);
GensIII1 := [];

CIII2 := HyperellipticCurve(6*x^5 + 14 * 5);
JIII2 := Jacobian(CIII2);
GensIII2 := [JIII2![x^2 - 8*x + 121/16, -435/8*x + 2191/32], JIII2![x - 3/2, -43/4], JIII2![x + 1, -8]];


CIII2 := HyperellipticCurve(6*x^5 + 14 * 5);
JIII2 := Jacobian(CIII2);
GensIII2 := [JIII2![x^2 - 8*x + 121/16, -435/8*x + 2191/32], JIII2![x - 3/2, -43/4], JIII2![x + 1, -8]];


// The case IV

CIV := HyperellipticCurve(30*x^5 + 70);
JIV := Jacobian(CIV);
GensIV := [JIV![x - 1, -10]];
