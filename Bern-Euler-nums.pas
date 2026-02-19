	{$MODE OBJFPC}{$H+}
	{$MINFPCONSTPREC 64}
	{$MODESWITCH ADVANCEDRECORDS}
	///{$CALLING pascal} 
	
	program test_unitDecfloat;
	
	uses sysutils, unitDecfloat;
	
	const N = 10;
		  BN = 41;
	var
		Xi: array [0..N] of decfloat;
		Wi: array [0..N] of decfloat;
		c: array [0..BN] of decfloat;
		b: array [0..BN] of decfloat;
		e: array [0..BN] of decfloat;
		y, z, re, one, two, tmp:decfloat;
		s:ansiString;
		ex, i, j, k:int32;

	//snippet borrowed from Rosetta Code https://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature#Pascal
	function F(X:decfloat) : decfloat;
	begin
	  F := Exp(X);
	end;
	 
	function LegInt(A,B:decfloat) : decfloat;
	var I      : integer;
		C1, C2: decfloat;
	begin
	  C1 := (B-A)/2;
	  C2 := (B+A)/2;
	  Result := 0;
	  For I := 1 to N do
		Result := Result + Wi[I] * F(C1*XI[I] + C2);
	  Result := C1 * Result;
	  LegInt := Result;
	end;
	// end snippet
	
	function fact(k : uint32) : decfloat;
	var
		f:decfloat;
		i:int32;
	begin
		f:=1;
		if k<2 then exit(f);
		for i:=2 to k do
			f:= f*i;
		fact:=f;
	end;
	
	begin
	gauss_leg_rule(N, Xi, Wi);
	writeln('   Gauss-Legendre ',N,' Quuadrature rule');
	writeln('                    Xi                                             Wi');
	for i:=1 to N do
	begin
		s:='   '+Xi[i].toString(42);
		while length(s)<50 do s:=s+' ';
		s:=s+Wi[i].toString(42);
		writeln(s);
	end;
	z:=LegInt(decfloat(-3), decfloat(3));
	y:=Exp(decfloat(3))-Exp(decfloat(-3));
	re:=(z-y)/y;
	ex:=abs(getExponent(re));
	if ex<4 then ex:=4;
	if ex>42 then ex:=42;
	Writeln;
	Writeln('Integrating Exp(x) over [-3, 3]: ',z.toString(ex));
	Writeln('Actual value: ',y.toString(ex));
	Writeln('Relative error: ',re.toString(10));
// =======================================================================================
	Writeln;
	one := 1;
	For i := 0 To BN do
	begin
		c[i] := one / (i+1);
		For j := i downto 1 do
		begin
			tmp := c[j - 1]-c[j];
			c[j - 1] := tmp * j;
		end;
		b[i] :=c[0];
		e[i] := 0;
	end;

	b[1]:=-b[1];
	For i := 0 To BN do
	begin
		s:=trim(inttoStr(i));
		while length(s)<3 do s:=' '+s;
		if i>2 then
		begin
			if not odd(i) then
				writeln(' Bernoulli(',s,') = ', b[i].toString(42));
		end
		else
			writeln(' Bernoulli(',s,') = ', b[i].toString(42));
	end;

	writeln;
	writeln( ' Euler numbers via Bernoulli numbers');
	writeln( ' https://math.stackexchange.com/questions/4489269/euler-numbers-and-bernoulli-numbers');

	writeln( '         n');
	writeln( '       -----');
	writeln( '        \       n!                          B(k+1)');
	writeln( ' E(n) =  )  --------- (2^(k+1) - 2^(2*k+2) -------');
	writeln( '        /    k!(b-k)!                        k+1');
	writeln( '       -----');
	writeln( '        k=0');
	writeln;

	two:=2;
	for i := 0 to BN-1 do
	begin
		for k:=0 to i do
		begin
			//k1:=1+k;
			//e(i)=e(i)+((two^k1-two^(2+2*k))*f1*b(k1))/((k1)*fact(k)*fact(i-k))
			e[i]:=e[i]+((fpipow(two,1+k)-fpipow(two, 2+2*k))*b[1+k]*fact(i))/((1+k)*fact(k)*fact(i-k));
		end;
	end;

	writeln(' Euler(  0) = ', e[0].toString(42));
	For i := 0 To BN do
	begin
		s:=trim(inttoStr(i));
		while length(s)<3 do s:=' '+s;
		if abs(e[i])>1 then
			writeln(' Euler(',s,') = ', e[i].toString(42));
	end;

end.
