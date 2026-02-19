	{$MODE OBJFPC}{$H+}
	{$MINFPCONSTPREC 64}
	{$MODESWITCH ADVANCEDRECORDS}
	///{$CALLING pascal} 
	
	program test_unitDecfloat;
	
	uses sysutils, windows, math, unitDecfloat;

	type
		TMatrix = array of array of DecFloat;
		TVector = array of DecFloat;

	{----------------------------------------------------------
	 Jacobi eigenvalue algorithm for symmetric matrix
	----------------------------------------------------------}
	procedure JacobiEigen(var A: TMatrix; var V: TMatrix; n: Integer);
	var
		i, j, p, q, iter: Integer;
		t, c, s, tau, temp, maxval, one: DecFloat;
	begin
		SetLength(V, n, n);
		
		one:=1;
		{ Initialize eigenvector matrix to identity }
		for i := 0 to n-1 do
			for j := 0 to n-1 do
				if i = j then
					V[i,j] := 1
				else
					V[i,j] := 0;
		
		for iter := 1 to 100*n*n do
		begin
			{ Find largest off-diagonal element }
			maxval := 0;
			p := 0; q := 1;
			
			for i := 0 to n-1 do
				for j := i+1 to n-1 do
					if Abs(A[i,j]) > maxval then
					begin
						maxval := Abs(A[i,j]);
						p := i;
						q := j;
					end;
			
			if maxval < 1e-14 then
				Break;
			
			tau := (A[q,q] - A[p,p]) / (2 * A[p,q]);
			
			if tau >= 0 then
				t := one / (tau + Sqrt(one + tau*tau))
			else
				t := -one / (-tau + Sqrt(one + tau*tau));
			
			c := one / Sqrt(one + t*t);
			s := t * c;
			//writeln(iter,'  ',t.toString);
			//writeln(iter,'  ',c.toString);
			{ Rotate matrix }
			for j := 0 to n-1 do
				if (j <> p) and (j <> q) then
				begin
					temp := A[p,j];
					A[p,j] := c*temp - s*A[q,j];
					A[j,p] := A[p,j];
					A[q,j] := s*temp + c*A[q,j];
					A[j,q] := A[q,j];
				end;
			
			temp := c*c*A[p,p] - 2*s*c*A[p,q] + s*s*A[q,q];
			A[q,q] := s*s*A[p,p] + 2*s*c*A[p,q] + c*c*A[q,q];
			A[p,p] := temp;
			A[p,q] := 0;
			A[q,p] := 0;
			
			{ Rotate eigenvectors }
			for j := 0 to n-1 do
			begin
				temp := V[j,p];
				V[j,p] := c*temp - s*V[j,q];
				V[j,q] := s*temp + c*V[j,q];
			end;
		end;
	end;

	{----------------------------------------------------------
	 Golub-Welsch Gauss-Hermite
	----------------------------------------------------------}
	procedure GaussHermiteGW(n: Integer);
	var
		i, j: Integer;
		A: TMatrix;
		V: TMatrix;
		nodes: TVector;
		weights: TVector;
		beta, sm: DecFloat;
	begin
		SetLength(A, n, n);
		
		{ Build Jacobi matrix }
		for i := 0 to n-1 do
			A[i,i] := 0;
		
		for i := 1 to n-1 do
		begin
			beta := i;
			beta := Sqrt(beta / 2);
			A[i-1,i] := beta;
			A[i,i-1] := beta;
		end;
		
		{ Compute eigenvalues and eigenvectors }
		JacobiEigen(A, V, n);
		
		{ Extract nodes and weights }
		SetLength(nodes, n);
		SetLength(weights, n);
		
		for i := 0 to n-1 do
		begin
			nodes[i] := A[i,i];
			weights[i] := V[0,i] * V[0,i] * Sqrt(fpPi);
		end;
		
		{ Sort (eigenvalues may not be ordered) }
		for i := 0 to n-2 do
			for j := i+1 to n-1 do
				if nodes[j] < nodes[i] then
				begin
					// Swap nodes
					beta := nodes[i];
					nodes[i] := nodes[j];
					nodes[j] := beta;
					
					// Swap weights
					beta := weights[i];
					weights[i] := weights[j];
					weights[j] := beta;
				end;
		
		{ Print }
		WriteLn('Gauss-Hermite Quadrature (Golub-Welsch)');
		WriteLn('n = ', n);
		WriteLn;
		WriteLn('   i               x_i                                          w_i');
		WriteLn('--------------------------------------------------------------------------------------');
		sm:=0;
		for i := 0 to n-1 do
		begin
			Write(Format('  %2d   ', [i+1]));
			writeln(nodes[i].toString(30),'       ', weights[i].toString(30));
			sm:=sm+weights[i];
		end;
		writeln;
		writeln('the sum of the weights = Sqrt(Pi)');
		writeln(sm.toString);
		
	end;

	{----------------------------------------------------------
	 MAIN
	----------------------------------------------------------}
	var
		n: Integer;
	begin
		Write('Enter order n: ');
		ReadLn(n);
		WriteLn;
		
		GaussHermiteGW(n);

		WriteLn;
		WriteLn('Done.');
		
		{ Pause to see output (optional) }
		if IsConsole then
		begin
			WriteLn;
			Write('Press Enter to exit...');
			ReadLn;
		end;
	end.

end.

