// Special cases for n = 5

// Curve CI5


function Hbound( : S := [2])
	/*
	INPUT:
		- S : a list of primes

	OUTPUT:
		An upper bound of the h(x) for an S-integral point P = (x,y)
	*/

	// We consider the equivelant curve y^2 = x^5 + 35 * 3^4 * 2^5
	R<x> := PolynomialRing(Rationals());
	K<rk> := NumberField(x^5 + 2^5*3^4*5*7);
	f := x^5 + 2^5*3^4*5*7;
	//theta := rk;
	//f := X^4 - theta*X^3 + theta^2*X^2 - theta^3*X + theta^4;
	C := HyperellipticCurve(x^5 + 35 * 3^4 * 2^5);
	
	J := Jacobian(C);
	GensJ := [J![x^2 - 48*x + 1089/4, 1305/4*x - 19719/8], J![x - 9, -387], J![x + 6, -288]];

	Sel,SeltoA := TwoSelmerGroup(J);
	
	H := [];
	for s in Sel do
		di := SeltoA(s);
		if #Coefficients(di) ne 1 then
			k := &+[x^(i-1) * Coefficients(di)[i] : i in [1..5]];
		else
			k := R!1;
		end if;
		//print "di",di,d;
		//print "k",k;
		hi := heightxBoundx(f,k,S);
		Append(~H,hi);
		//printf "hi = %o\n",hi;
	end for;

	Hmax := Max(H);	
	print "max(H)",Hmax;	

	return Hmax;
end function;



function apply_MW_sieve(J,mwBasis,Hbound : primebound := 10^4)
	/*
	INPUT:
		- J : the Jacobian of C
		- mwBasis : generators of J(Q)
		- Hbound : an upper bound of h(x) for a S-integral point P = (x,y)
		- primebound : an upper bound for the primes we use in the Mordell-Weil sieve
	*/

	C := Curve(J);
	
	// Both B and W are given by the MordellWeillSieve function
	B := 4449329780614748206472972686179940652515754483274306796568214048000;
	W := [P - C![1,0,0] : P in Points(C : Bound := 1000)];
	
	mu1 := HeightConstant(J: Factor:=true, Effort:=2);
	mu1 := -mu1;
	mu2 := Max([Sqrt(Height(i)):i in W]);
	mu3:= Min([Sqrt(l[1]) : l in Eigenvalues(HeightPairingMatrix(mwBasis))]);
	
	print "mu1 is ", mu1;
	print "mu2 is ", mu2;
	print "mu3 is ", mu3; 


	F := FreeAbelianGroup(3);
	L := sub<F | [B*F.1,B*F.2,B*F.3]>;
	
	L,B,iL,mL,successLevels,bprimes,lprime := mwSievePrimesEnhance(F,L,W,J,mwBasis,11,primebound,B,[2,3,5,7]);
	

	assert(mu3*mL-mu2 gt 0);	
	lowerBound := (mu3*mL-mu2)^2 + mu1;
	print "lowerBound",lowerBound;
	
	return true;

end function;



// Ruiz code


//////////////////////////////////////////////////////////////
// Procedures in the context of                             //
// S-integral points on hyperelliptic curves                //
// 10 March 2010                                            //
//////////////////////////////////////////////////////////////
// Computation of the constants in Sections 4-8             //
//////////////////////////////////////////////////////////////

// The folowing function
// returns c_1(s,d),...,c_5(s,d),
// the constants from Lemmas 4.1, 4.2, 4.3.
cConstants:=function(s,d);
        c1:=Factorial(s-1)^2/(2^(s-2)*d^(s-1))+0.0;
        c2:=29*Exp(1)*Sqrt(s-2)*d^(s-1)*Max(Log(d),1)*c1;
        if d eq 1 then;
            c3:=Factorial(s-1)^2/2^(s-1)*2/Log(2)+0.0;
        else
            c3:=Factorial(s-1)^2/2^(s-1)*Log(3*d)^3+0.0;
        end if;
        c4:=d*Pi(RealField())^(s-2)*c2;
        c5:=2*d*c3;
        return c1,c2,c3,c4,c5;
end function;

// The following function
// returns c6(r,d),
// the constant from Lemma 4.4
cConstants1:=function(r,d);
	if r eq 0 then;
	   c6:=0.0;
	else
	   if r eq 1 then;
	      c6:=1/d+0.0;
	   else
	      c6:=29*Exp(1)*Factorial(r)*r*Sqrt(r-1)*Log(d)+0.0;
	   end if;
	 end if;
	 return c6;
end function;

// The following function
// returns c7,c8,c9,
// the constants from Lemmas 6.1 and 6.2
cConstants2:=function(n,d);
        c7:=Min(1.451*(30*Sqrt(2))^(n+4)*(n+1)^(5.5),Pi(RealField())*2^(6.5*n+27))*d^2*Log(Exp(1)*d);
	c8:=(16*Exp(1)*d)^(2*(n+1))*n^(1.5)*Log(2*n*d)*Log(2*d);
	c9:=(2*d)^(2*n+1)*Log(2*d)*Log(3*d)^3;
	return c7,c8,c9;
end function;

// The following function returns an
// upper bound for the regulator times the class number
// of a number field K and upper bound for the $S$-regulator
// coming from Lemmas 5.1 and 5.2
// Input: 
// u the number of real places of K
// v the number of complex places
// L an upper bound for the discriminant.
// w the number of roots of unity in K
// P largest prime below the places in $S$
// t number of finite places in S
regBound:=function(u,v,L,w,P,t);
	d:=u+2*v;
	a:=2^(-v)*Sqrt(L/Pi(RealField())^d);
	f:=2^(-u)*w*a^2*2^(d+1)+0.0;
        for j in [1..999] do
                i:=2-j/1000.0;
                f1:=2^(-u)*w*a^i*Gamma(i/2)^u*Gamma(i)^v*i^(d+1)*(i-1)^(1-d);
                f:=Min(f,f1);
        end for;
	if t eq 0 then
            return f,f;
	else
	    return f,f*(d*Max(Log(P),1))^t;
        end if;
end function;

// The following function returns an
// upper bound for the class number field K
// coming from Lemma 5.3
// Input:
// L an upper bound for the discriminant of K
// u the number of real places of K
// v the number of complex places of K
classNumberBound:=function(L,u,v);
	d:=u+2*v;
	a:=(2/Pi(RealField()))^v*Sqrt(L);
	return a*(d-1+Log(a))^(d-1)/Factorial(d-1);
end function;

// This function gives a bound on $x$ such that
// x-\alpha = \kappa \xi^2
// The bound is the one given in Theorem 8.2.
// input p, k,
// p the monic polynomial defining \alpha
// k a polynomial s.t. k(\alpha)= \kappa
// S a finite set of rational primes
heightxBoundx:=function(p,k,S);
        d:=Degree(p);
        Qa1<alpha1>:=NumberField(p); // QQ(alpha1).
	Qz<z>:=PolynomialRing(Qa1);  // We'll factor p over Qa1(z).
	pz:=Qz!p;
	p1:=Factorisation(pz)[2][1];  // This is the non linear factor
	Qa2<alpha2>:=NumberField(p1); // such that one of its roots is alpha2.
	Qa12<a>:=AbsoluteField(Qa2);  // Qa2 as extension over QQ.
        if Degree(Qa12) ne d*(d-1) then
                error("program written on assumption that Galois action is 2-transitive"); // this makes K_1, K_2, K_3 in Thm 2 isomorphic 
        end if;
        OO:=MaximalOrder(Qa12);
        OOy<y>:=PolynomialRing(OO);
        kappa1:=OO!(Qa12!Evaluate(k,alpha1));
        kappa2:=OO!(Qa12!Evaluate(k,alpha2));
        kappa12:=kappa1*kappa2;
        if #Factorisation(y^2-kappa12) eq 1 then
                subOrderK:=ext<OO | y^2-kappa12>;
                subOrderK:=AbsoluteOrder(subOrderK);
                D:=Discriminant(subOrderK);
                for p in PrimeDivisors(D) do
                        subOrderK:=pMaximalOrder(subOrderK,p);
                        // remarkably this gives the maximal order much
                        // quicker than MaximalOrder();
                end for;
                OOK:=subOrderK; // maximal order
                dL:=4*d*(d-1)*(d-2);  //upper bound for the degree of L
        else
                OOK:=OO;
                dL:=d*(d-1)*(d-2);
        end if;
       K:=NumberField(OOK);
       dK:=Degree(K);
       L:=AbsoluteValue(Discriminant(OOK));
        u,v:=Signature(K);
        s:=u+v;
	SumLogsPrimes:=0.0;
        for p in S do
              s:=s + #Decomposition(OOK,p); //Compute the number of
                                          // places over the primes in S.
	      SumLogsPrimes:=SumLogsPrimes+Log(p); //This will appear in Hs
        end for;
        r:=u+v-1;//Unit rank
        if u gt 0 then  // w will be the number of roots of unity
                        //
                w:=2;   // MAGMA doesn't seem to check if a field has
                        // real embeddings before it computes the 
                        // roots of unity!!
                        // and it is very slow for fields with real embd!!
        else
                w:=#TorsionUnitGroup(OOK);
        end if;
      P:=Max(Append(S,1));  // Largest prime in S (1 if S is empty)
      cs1,cs2,cs3,cs4,cs5:=cConstants(s,dK); // $c^*_j(s_i,d_i)$
      t:=s-(r+1);           // Number of finite places over primes in S
      cs6:=cConstants1(r,dK);
      c7,c8,c9:=cConstants2(2*s-1,dL);
      R,RS:=regBound(u,v,L,w,P,t); //upper bound for regs and S-regs of K_i
      N:=(Norm((Qa12!kappa1)*((Qa12!alpha1)-(Qa12!alpha2))))^2;
      hkappa:=AbsoluteLogarithmicHeight(kappa1);
      halpha:=AbsoluteLogarithmicHeight(alpha1);
      h:=classNumberBound(L,u,v);
      Hs:=Max(cs6*R+Log(N)/dK+h*SumLogsPrimes+hkappa,1);
      Hs:=Max(Hs,Pi(RealField())/dL);
      cs10:=2*Hs+2*Hs*dL*(#S+1)*(1+2*cs4^2)*c7*RS^2*Log(Sqrt(2)*Exp(1)*Max(cs2*RS,(2*s-1)*Pi(RealField())/Sqrt(2)));
      cs11:=4*dL*(#S+1)*Hs*cs4^2*c7*RS^2;
      cs12:=2*Hs+2*Hs*(dL*(#S+1))+cs11*Log(Max(cs5,1)/(2*Sqrt(2)*dL*Hs));
      cs13:=Log(2)+2*Hs+4*(2*s-2)*Hs*cs1^2*cs2*c9*RS^3;
      cs14:=2*Hs*dL^(2*s-2)*P^dL*cs1^2*c8*RS^2/(Log(2)*Max(1,Log(P^dL)));
      cs15:=2*Hs+2*Hs*dL*(#S+1)+cs14*Log(Max(cs5,1)*Exp(2*s*(12*s-1))*dL^(3*(2*s-1))*Log(2*dL)*P^(2*s*dL)/(Hs*c9));
      B:=20*Log(2)+13*hkappa+19*halpha+Hs;
      M:={cs10/2,cs12+cs11*Log(cs11),cs13/2,cs15+cs14*Log(cs14)};
//      print M;
      Minfinite:={cs10/2,cs12+cs11*Log(cs11)};
      if #S eq 0 then;
            B:=B+8*Max(Minfinite);
      else
	    B:=B+8*Max(M);
      end if;
      print "Bound for the S-reg", RS;
      print "S-unit rank", s;
      print "bound for this kappa", B;
      print "====================";
      return B;
end function;




// Mordell-Weil sieve



// This function attempts to replace $L$ be a sublattice
// as in [11, Lemma 11.1]

// Here J is the Jacobian
// W is the image on J of the known rational points on C
// mwBasis is the basis for the free part of J(Q)
// q is a prime of good reduction
// F:=FreeAbelianGroup(r); where $r$ is the mw rank
// L is some subgroup of F
// B is some integer such that L \subseteq B \Z^r

// returns tf, n, L
// where tf is true if all of I--IV succeed and false otherwise
// n is the success level, 1--4 if it fails on I--IV and 5 if it succeeds on all
// L is the old lattice if failure and the new lattice if success.
mwSievePrimeEnhance:=function(F,L,W,J,mwBasis,q,B);
	r:=#mwBasis;
	Jq:=BaseChange(J,GF(q));

	A,psi:=AbelianGroup(Jq); //computing the group structure of J mod q;
	ordA:=#A;

	if GCD(B,ordA) lt (ordA^0.6) then
		//print "Criterion I fails";
		return false, 1, L;
	end if;
	
	mwRed:=[(Jq!D)@@psi : D in mwBasis];

	phi:=hom<F->A | mwRed>;
	phires:=hom<L->A| [(F!(L.i))@phi : i in [1..r]]>;
	
	Lnew:=Kernel(phires);

	Q,theta:=quo<L|Lnew>;

	if #Q eq 1 then
		//print "Criterion II fails";
		return false, 2, L;
	end if;

	reps:=[F!(q@@theta) : q in Q | q ne Q!0];  // non-zero coset reps for L/Lnew
	notests:=(#W)*(#reps);  // Approximate number of tests to do
	if notests gt 2*q then
		//print "Criterion III fails";
		return false, 3, L;
	end if;

	for w in W do
		wim:=(Jq!w);
		for l in reps do
			P:=wim+(l@phi)@psi;
			if Degree(P[1]) le 1 then
				//print "Criterion IV fails";
				return false, 4, L;
			end if;
		end for;
	end for;
	Lnew:=sub<F | [F!(L!Lnew.i) : i in [1..r]]>; //identify Lnew 
						// as subgroup of F
	print "Criteria I--IV succeeded for prime", q;
	print "the new index is", Index(F,Lnew)+0.0; // printing index as 
							// a real number
							// o/w it is too long
	return true, 5, Lnew;
end function;

// repeatedly apply the above using all primes from intialprime up to the n-th prime
// where n = primebound,
// badprimes is the set of primes that are of bad reduction, or primes we don't
// check for they have been already used.
// At the end we compute a new constant newB such that the L is contained in newB Z^r
// so that we can repeat the process but using primes up to the n-th prime that did
// not give information in the previous application.
mwSievePrimesEnhance:=function(F,L,W,J,mwBasis,initialprime,primebound,B,badprimes);
        r:=#mwBasis; //rank of J
        amb:=StandardLattice(r);
        successLevels:=[0,0,0,0,0]; // this records the number of primes
                                    // that fail at criteria 1,..,4
                                    // or succeed in all 
                                    // first entry is number of those that fail at criterion 1
                                    // last is number of those that succeed in all criteria
        Llattice:=amb;
        prime:=initialprime;  // 'prime' is the current prime we are working with
      for q:=1 to primebound do;
            if (prime notin badprimes) then; 
                // print "using the prime", prime, "to enhance the information";
                tf, level, L1:=mwSievePrimeEnhance(F,L,W,J,mwBasis,prime,B);
                successLevels[level]:=successLevels[level]+1;
                if tf then
                        badprimes:=Append(badprimes,prime); // the prime will not give us 
                                                            // new information in the next pass
                        L:=L1;
                        Llattice:=sub<amb | [Eltseq(F!(L.i)) : i in [1..r]]>;
                        mL:=Length(ShortestVectors(Llattice)[1]);
                        print "the shortest vector of the new lattice L has length", mL;
//                   print "success levels so far", successLevels;
                end if;
          end if;
            prime:=NextPrime(prime);
        end for;
      newB:=0;
      for g in Generators(L) do;
                newB:=GCD(GCD(Eltseq(F!g)),newB); 
      end for;
 // We will have now j(C) is contained in W + newB*J
return L,newB,Index(F,L),Length(ShortestVectors(Llattice)[1]),successLevels,badprimes,prime;
end function;


