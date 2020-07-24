// We compute classical newforms using the Hilbert newforms package


function Hilbert_newforms(N)
	/*
	INPUT:
		- ``N`` : an integer

	OUTPUT:
		The space of classical newforms of level N as Hilbert newforms
	*/

	QQ := RationalsAsNumberField();
	ZZ := Integers(QQ);
	M := HilbertCuspForms(QQ, N*ZZ);
	p := Min([p : p in PrimeFactors(N) | Valuation(N,p) eq 1]);
	A := QuaternionAlgebra(p*ZZ, InfinitePlaces(QQ) : Optimized);
	NewFs := NewSubspace(M : QuaternionOrder:=MaximalOrder(A));

	return  NewformDecomposition(NewFs);

end function;



function Dimension_Hilbert(NewFs)

	dim := 0;

	for fnew in NewFs do
		dim +:= Dimension(fnew);
	end for;



	n := #NewFs;
	maxd := Dimension(NewFs[n]);

	degs := [0 : i in [1..maxd]];

	for fnew in NewFs do
		degfnew := Dimension(fnew);
		degs[degfnew] +:= degfnew;
	end for;


	return dim,#NewFs,degs;
end function;



