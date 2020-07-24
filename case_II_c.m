/* We work on the case II subcase c*/


function case_II_c( : bound := 20, NewE2s := [], NewF2s := [])
	/*
	INPUT:
		- ``NewE2s`` : the space of newforms of level 2^8 * 5 * 7 as HIlbert newforms.
		- ``NewF2s`` : the space of newforms of level 2^3 * 3 * 5^2 * 7 as HIlbert newforms.

	OUTPUT:
		True, if the elimination step does not fail together with a list of small exponents.
	*/


	if #NewE2s eq 0 then
		NewE2s := Hilbert_newforms(2^8 * 5 * 7);
	end if;


	if #NewF2s eq 0 then
		NewF2s := Hilbert_newforms(2^3 * 3 * 5^2 * 7);
	end if;

	primes_q := [q : q in PrimesUpTo(bound) | q ge 11];
	ZZE2 := RingOfIntegers(BaseField(NewE2s[1]));
	ZZF2 := RingOfIntegers(BaseField(NewF2s[1]));

	
	i := 1;
	for E2new in NewE2s do
		for F2new in NewF2s do
			print "pair : ",i;
			i +:= 1;
			Ppair := [* *];
			for q in primes_q do
				aqE2new := HeckeEigenvalue(Eigenform(E2new),q*ZZE2);
				aqF2new := HeckeEigenvalue(Eigenform(F2new),q*ZZF2);


				Zq := FiniteField(q);
				Pq := q;
				for y1 in [0..q-1] do
					for a in [0..q-2] do
						for b in [0..q-2] do
							z1 := y1^2;

							// F22 Frey curve
							AF2 :=   - (3 * z1* 5^(2*b) + 40);
							BF2 := 280;

							// The discriminant of F2
							DF2 := BF2^2 * (-AF2^2 + 4*BF2);


							//Frey curve for E2
							AE2 := 4 * (z1* 5^(2*b) + 4);
							BE2 :=  14* 5^(4*b-1)*z1^2;

							// The discriminant of E2
							DE2 := BE2^2 * (-AE2^2 + 4*BE2);

							DF2 := Zq!DF2;
							DE2 := Zq!DE2;						

							if DF2 ne 0 and DE2 ne 0 then
								aqF2 := TraceOfFrobenius(EllipticCurve([0, AF2,0,BF2,0]),q);
								aqE2 := TraceOfFrobenius(EllipticCurve([0, AE2,0,BE2,0]),q);

								dE2 := aqE2new - aqE2;
								dF2 := aqF2new - aqF2;
						
							elif DF2 eq 0 and DE2 ne 0 then
								aqE2 := TraceOfFrobenius(EllipticCurve([0, AE2,0,BE2,0]),q);

								dE2 := aqE2new - aqE2;
								dF2 := (q+1)^2 - aqF2new^2;
							
							elif DF2 ne 0 and DE2 eq 0 then
								aqF2 := TraceOfFrobenius(EllipticCurve([0, AF2,0,BF2,0]),q);

								dE2 := (q+1)^2 - aqE2new^2;
								dF2 := aqF2new - aqF2;
								
							else
								dE2 := (q+1)^2 - aqE2new^2;
								dF2 := (q+1)^2 - aqF2new^2;
							end if;

							
							dE2 := Integers()!Norm(dE2);
							dF2 :=Integers()! Norm(dF2);
							
							if dE2 eq 0 and dF2 eq 0 then
								Pq *:= 0;
							elif dE2 ne 0 and dF2 eq 0 then
								Pq *:= dE2;
							elif dE2 eq 0 and dF2 ne 0 then
								Pq *:= dF2;
							else
								Pq *:= GCD(dE2,dF2);
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
						elimination_pair := case_II_c( : bound := more_qs, NewE2s := [*E2new*], NewF2s := [*F2new*]);
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
