// This function runs all case



load "/home/angelos/Magma/5th_powers/case_I_a.m";
load "/home/angelos/Magma/5th_powers/case_I_b.m";
load "/home/angelos/Magma/5th_powers/case_I_c.m";
load "/home/angelos/Magma/5th_powers/case_II_a.m";
load "/home/angelos/Magma/5th_powers/case_II_b.m";
load "/home/angelos/Magma/5th_powers/case_II_c.m";
load "/home/angelos/Magma/5th_powers/case_III_a.m";
load "/home/angelos/Magma/5th_powers/case_III_b.m";
load "/home/angelos/Magma/5th_powers/case_IV_a.m";
load "/home/angelos/Magma/5th_powers/case_IV_b.m";



// We compute spaces of Newforsm

load "/home/angelos/Magma/5th_powers/Hilbert_newforms.m";

Newf1011 := Hilbert_newforms(2 * 5 * 7);
Newf1111 := Hilbert_newforms(2 * 3 * 5 * 7);
Newf1021 := Hilbert_newforms(2 * 5^2 * 7);
Newf3111 := Hilbert_newforms(2^3 * 3 * 5 * 7);
Newf1121 := Hilbert_newforms(2 * 3 * 5^2 * 7);
Newf5011 := Hilbert_newforms(2^5 * 5 * 7);
Newf3121 := Hilbert_newforms(2^3 * 3 * 5^2 * 7);
Newf5021 := Hilbert_newforms(2^5 * 5^2 * 7);
Newf8011 := Hilbert_newforms(2^8 * 5 * 7);


// Most time expensive

time Newf7111 := Hilbert_newforms(2^7 * 3 * 5 * 7);
time Newf8111 := Hilbert_newforms(2^8 * 3 * 5 * 7);
time Newf7121 := Hilbert_newforms(2^7 * 3 * 5^2 * 7);
time Newf8021 := Hilbert_newforms(2^8 * 5^2 * 7);
time Newf8121 := Hilbert_newforms(2^8 * 3 * 5^2 * 7);
