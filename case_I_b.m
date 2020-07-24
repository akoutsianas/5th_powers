/* We work on the case I subcase b*/


function case_I_b( : bound := 20, NewE1s := [], NewF1s := [])
	/*
	INPUT:
		- ``NewE1s`` : the space of newforms of level 2^8 * 5^2 * 7 as Hilbert newforms
		- ``NewF1s`` : the space of newforms of level 2^7 * 3 * 5 * 7 as Hilbert newforms.

	OUTPUT:
		True, if the elimination step does not fail and p < 7.
	*/


	if #NewE1s eq 0 then
		NewE1s := Hilbert_newforms(2^8 * 5^2 * 7);
	end if;


	if #NewF1s eq 0 then
		NewF1s := Hilbert_newforms(2^7 * 3 * 5 * 7);
	end if;

	primes_q := [q : q in PrimesUpTo(bound) | q ge 11];
	ZZE1 := RingOfIntegers(BaseField(NewE1s[1]));
	ZZF1 := RingOfIntegers(BaseField(NewF1s[1]));

	primes_q := [q : q in PrimesUpTo(bound) | q ge 11];


	exc_pairs := [];
	i := 1;
	for E1new in NewE1s do
		for F1new in NewF1s do
			
			print "pair",i;
			i +:= 1;
			Ppair := [* *];
			for q in primes_q do
				aqE1new := HeckeEigenvalue(Eigenform(E1new),q*ZZE1);
				aqF1new := HeckeEigenvalue(Eigenform(F1new),q*ZZF1);


				Zq := FiniteField(q);
				Pq := q;
				for y1 in [0..q-1] do
					for a in [0..q-2] do
						for b in [0..q-2] do
							z1 := y1^2;

							// F1 Frey curve
							AF1 := 2 * (3*z1 + 2 * 5^(2*b + 1));
							BF1 := 14 * 5^(4*b + 1);

							// The discriminant of F1
							DF1 :=  BF1^2 * (-AF1^2 + 4*BF1);


							//Frey curve for E1
							AE1 := 20 * (z1 + 5^(2*b));
							BE1 :=  70 * z1^2;

							// The discriminant of E1
							DE1 := BE1^2 * (-AE1^2 + 4*BE1);

							DF1 := Zq!DF1;
							DE1 := Zq!DE1;						

							if DF1 ne 0 and DE1 ne 0 then
								aqF1 := TraceOfFrobenius(EllipticCurve([0, AF1,0,BF1,0]),q);
								aqE1 := TraceOfFrobenius(EllipticCurve([0, AE1,0,BE1,0]),q);

								dE1 := aqE1new - aqE1;
								dF1 := aqF1new - aqF1;
						
							elif DF1 eq 0 and DE1 ne 0 then
								aqE1 := TraceOfFrobenius(EllipticCurve([0, AE1,0,BE1,0]),q);

								dE1 := aqE1new - aqE1;
								dF1 := (q+1)^2 - aqF1new^2;
							
							elif DF1 ne 0 and DE1 eq 0 then
								aqF1 := TraceOfFrobenius(EllipticCurve([0, AF1,0,BF1,0]),q);

								dE1 := (q+1)^2 - aqE1new^2;
								dF1 := aqF1new - aqF1;
								
							else
								dE1 := (q+1)^2 - aqE1new^2;
								dF1 := (q+1)^2 - aqF1new^2;
							end if;

							
							dE1 := Integers()!Norm(dE1);
							dF1 :=Integers()! Norm(dF1);
							
							if dE1 eq 0 and dF1 eq 0 then
								Pq *:= 0;
							elif dE1 ne 0 and dF1 eq 0 then
								Pq *:= dE1;
							elif dE1 eq 0 and dF1 ne 0 then
								Pq *:= dF1;
							else
								Pq *:= GCD(dE1,dF1);
							end if;

						end for;
					end for;
				end for;				
				
				Append(~Ppair,Pq);
			end for;

			if #[dif: dif in Ppair | dif ne 0] ne 0 then
				Pr := PrimeFactors(GCD([dif: dif in Ppair | dif ne 0]));
				
				if #Pr gt 0 and Max(Pr) gt 5 then
					
					more_qs := bound + 20;
					if more_qs le 100 then
						elimination_pair := case_I_b( : bound := more_qs, NewE1s := [*E1new*], NewF1s := [*F1new*]);
					else
						return false;
					end if;

					if not elimination_pair then
						print "We need more primes for a pair to get n<= 5";
						return false;
					end if;
				end if;
			else
				print "We have zeros for all q";
				return false;
			end if; 
		end for;
	end for;	


	return true;
end function;
