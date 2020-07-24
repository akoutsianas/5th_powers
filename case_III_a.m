/* We worgenvalue(Eigenform(F3new),q*ZZF3);k on the case III subcase a*/


function case_III_a( : bound := 20, NewE3s := [], NewF3s := [])
	/*
	INPUT:
		- ``NewE3s`` : the space of newforms of level 2 * 5^2 * 7 as Hilbert newforms.
		- ``NewF3s`` : the space of newforms of level 2^8 * 3 * 5 * 7 as Hilbert newforms.

	OUTPUT:
		True, if the elimination step does not fail and p < 7.
	*/


	if #NewE3s eq 0 then
		NewE3s := Hilbert_newforms(2 * 5^2 * 7);
	end if;


	if #NewF3s eq 0 then
		NewF3s := Hilbert_newforms(2^8 * 3 * 5 * 7);
	end if;

	primes_q := [q : q in PrimesUpTo(bound) | q ge 11];
	ZZE3 := RingOfIntegers(BaseField(NewE3s[1]));
	ZZF3 := RingOfIntegers(BaseField(NewF3s[1]));

	i := 1;
	for E3new in NewE3s do
		for F3new in NewF3s do
			print "pair : ",i;
			i +:= 1;
			Ppair := [* *];
			for q in primes_q do
				aqE3new := HeckeEigenvalue(Eigenform(E3new),q*ZZE3);
				aqF3new := HeckeEigenvalue(Eigenform(F3new),q*ZZF3);


				Zq := FiniteField(q);
				Pq := q;
				for y1 in [0..q-1] do
					for a in [0..q-2] do
						for b in [0..q-2] do
							z1 := y1^2;

							// F1 Frey curve
							AF3 :=  4*(3 * 2^(2*a - 1) * z1 + 5^(2*b + 1));
							BF3 := 14 * 5^(4*b+1);

							// The discriminant of F1
							DF3 := BF3^2 * (4*BF3 - AF3^2);


							//Frey curve for E1
							AE3 := (5 * (2^(2*a) * z1 + 5^(2*b)) - 1)/4;
							BE3 := 35 * 2^(4*a - 7) * z1^2;

							// The discriminant of E1
							DE3 := BE3^2 * (-16*AE3^2 - 8*AE3 + 64*BE3 - 1);

							DF3 := Zq!DF3;
							DE3 := Zq!DE3;						

							if DF3 ne 0 and DE3 ne 0 then
								aqF3 := TraceOfFrobenius(EllipticCurve([0, AF3,0,BF3,0]),q);
								aqE3 := TraceOfFrobenius(EllipticCurve([1, AE3,0,BE3,0]),q);

								dE3 := aqE3new - aqE3;
								dF3 := aqF3new - aqF3;
						
							elif DF3 eq 0 and DE3 ne 0 then
								aqE3 := TraceOfFrobenius(EllipticCurve([1, AE3,0,BE3,0]),q);

								dE3 := aqE3new - aqE3;
								dF3 := (q+1)^2 - aqF3new^2;
							
							elif DF3 ne 0 and DE3 eq 0 then
								aqF3 := TraceOfFrobenius(EllipticCurve([0, AF3,0,BF3,0]),q);

								dE3 := (q+1)^2 - aqE3new^2;
								dF3 := aqF3new - aqF3;
								
							else
								dE3 := (q+1)^2 - aqE3new^2;
								dF3 := (q+1)^2 - aqF3new^2;
							end if;

							
							dE3 := Integers()!Norm(dE3);
							dF3 :=Integers()! Norm(dF3);
							
							if dE3 eq 0 and dF3 eq 0 then
								Pq *:= 0;
							elif dE3 ne 0 and dF3 eq 0 then
								Pq *:= dE3;
							elif dE3 eq 0 and dF3 ne 0 then
								Pq *:= dF3;
							else
								Pq *:= GCD(dE3,dF3);
							end if;

						end for;
					end for;
				end for;				
				
				Append(~Ppair,Pq);
			end for;
			
			//print "Ppair",#[dif: dif in Ppair | dif ne 0];

			if #[dif: dif in Ppair | dif ne 0] ne 0 then
				Pr := PrimeFactors(GCD([dif: dif in Ppair | dif ne 0]));
				if #Pr gt 0 and Max(Pr) gt 5 then

					more_qs := bound + 20;
					if more_qs le 100 then
						elimination_pair := case_III_a( : bound := more_qs, NewE3s := [*E3new*], NewF3s := [*F3new*]);
					else
						return false;
					end if;

					if not elimination_pair then
						print "We need more primes for a pair to get n<= 5 and Pr = ",Pr;
						return true;
					end if;
				end if;
			else
				more_qs := bound + 20;
				if more_qs le 100 then
					elimination_pair := case_III_a( : bound := more_qs, NewE3s := [*E3new*], NewF3s := [*F3new*]);
				else
					return false;
				end if;

				if not elimination_pair then
					print "We have zeros for all q";
					return false;
				end if;
			end if;
		end for;
	end for;	


	return true;
end function;
