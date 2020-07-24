/* We work on the case IV subcase a*/


function case_IV_b( : bound := 20, NewE4s := [], NewF4s := [])
	/*
	INPUT:
		- ``NewE4s`` : the space of newforms of level 2^5 * 5 * 7 as Hilbert newforms
		- ``NewF4s`` : the space of newforms of level 2^8 * 3 * 5^2 * 7 as Hilbert newforms

	OUTPUT:
		True, if the elimination step does not fail and p < 7.
	*/


	if #NewE4s eq 0 then
		NewE4s := Hilbert_newforms(2^5 * 5 * 7);
	end if;


	if #NewF4s eq 0 then
		NewF4s := Hilbert_newforms(2^8 * 3 * 5^2 * 7);
	end if;

	primes_q := [q : q in PrimesUpTo(bound) | q ge 11];
	ZZE4 := RingOfIntegers(BaseField(NewE4s[1]));
	ZZF4 := RingOfIntegers(BaseField(NewF4s[1]));


	i := 1;
	for E4new in NewE4s do
		for F4new in NewF4s do
			print "pair i : ",i;
			i +:= 1;
			Ppair := [* *];
			for q in primes_q do
				aqE4new := HeckeEigenvalue(Eigenform(E4new),q*ZZE4);
				aqF4new := HeckeEigenvalue(Eigenform(F4new),q*ZZF4);

				Zq := FiniteField(q);
				Pq := q;
				for y1 in [0..q-1] do
					for a in [0..q-2] do
						for b in [0..q-2] do
							z1 := y1^2;

							// F1 Frey curve
							AF4 :=  20 * (6*5^(2*b-1) * z1 + 1);
							BF4 := 70;

							// The discriminant of F1
							DF4 := BF4^2 * (-AF4^2 + 4*BF4);


							//Frey curve for E1
							AE4 := 2^2 * 5^(2*b) * z1 + 1;
							BE4 := 14 * 5^(4*b-1) * z1^2;

							// The discriminant of E1
							DE4 :=  BE4^2 * (-AE4^2 + 4*BE4);

							DF4 := Zq!DF4;
							DE4 := Zq!DE4;			

							if DF4 ne 0 and DE4 ne 0 then
								aqF4 := TraceOfFrobenius(EllipticCurve([0, AF4,0,BF4,0]),q);
								aqE4 := TraceOfFrobenius(EllipticCurve([0, AE4,0,BE4,0]),q);

								dE4 := aqE4new - aqE4;
								dF4 := aqF4new - aqF4;
						
							elif DF4 eq 0 and DE4 ne 0 then
								aqE4 := TraceOfFrobenius(EllipticCurve([0, AE4,0,BE4,0]),q);

								dE4 := aqE4new - aqE4;
								dF4 := (q+1)^2 - aqF4new^2;
							
							elif DF4 ne 0 and DE4 eq 0 then
								aqF4 := TraceOfFrobenius(EllipticCurve([0, AF4,0,BF4,0]),q);

								dE4 := (q+1)^2 - aqE4new^2;
								dF4 := aqF4new - aqF4;
								
							else
								dE4 := (q+1)^2 - aqE4new^2;
								dF4 := (q+1)^2 - aqF4new^2;
							end if;

							
							dE4 := Integers()!Norm(dE4);
							dF4 := Integers()!Norm(dF4);
							
						
							if dE4 eq 0 and dF4 eq 0 then
								Pq *:= 0;
							elif dE4 ne 0 and dF4 eq 0 then
								Pq *:= dE4;
							elif dE4 eq 0 and dF4 ne 0 then
								Pq *:= dF4;
							else
								Pq *:= GCD(dE4,dF4);
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
					if more_qs le 50 then
						elimination_pair := case_IV_b( : bound := more_qs, NewE4s := [*E4new*], NewF4s := [*F4new*]);
					else
						return false;
					end if;

					if not elimination_pair then
						print "We need more primes for a pair to get n<= 5";
						//return false;
					end if;
				end if;
			else
				more_qs := bound + 20;
				if more_qs le 100 then
					elimination_pair := case_IV_b( : bound := more_qs, NewE4s := [*E4new*], NewF4s := [*F4new*]);
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
