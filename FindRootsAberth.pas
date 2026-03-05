program AberthRoots;

{$mode objfpc}{$h+}
{$minfpconstprec 64}
{$modeswitch advancedrecords}

uses
	SysUtils, Math, unitDecfloat;

function C(r, i: DecFloat): DecFloatC;
begin
	Result.re := r;
	Result.im := i;
end;

Function AbsC(z : DecFloatC) : DecFloat;
begin
    Result := Sqrt(z.re*z.re + z.im*z.im);
End;

function PolyEval(const coef: array of DecFloat; const z: DecFloatC): DecFloatC;
var
	i: Integer;
	sum: DecFloatC;
begin
	sum := C(coef[High(coef)],0);
	for i := High(coef)-1 downto 0 do
		sum := sum * z + C(coef[i],0);
	Result := sum;
end;

function PolyDeriv(const coef: array of DecFloat; const z: DecFloatC): DecFloatC;
var
	i: Integer;
	sum: DecFloatC;
begin
	if High(coef) <= 0 then
		Exit(C(0,0));

	sum := C(High(coef)*coef[High(coef)],0);

	for i := High(coef)-1 downto 1 do
		sum := sum * z + C(i*coef[i],0);

	Result := sum;
end;

procedure AberthInitialGuesses(var roots: array of DecFloatC; n: Integer);
var
	k: Integer;
	r, theta: DecFloat;
begin
	Randomize;
	r := 1.2 + Random * 0.8;

	for k := 0 to n-1 do
	begin
		theta := 2 * Pi * k / n + Random * 0.4;
		roots[k].re := r * Cos(theta);
		roots[k].im := r * Sin(theta);
	end;
end;

procedure FindRootsAberth(const coef: array of DecFloat; var roots_out: array of DecFloatC);
var
	n, k, j, iter, maxiter, converged: Integer;
	roots: array of DecFloatC;
	zk, pz, dpz, correction_inv, sigma, diff, w: DecFloatC;
begin
	n := High(coef);
	if n < 1 then Exit;

	SetLength(roots, n);

	AberthInitialGuesses(roots, n);

	maxiter := 120;
	iter := 0;
	converged := 0;

	while (converged < n) and (iter < maxiter) do
	begin
		Inc(iter);
		converged := 0;

		for k := 0 to n-1 do
		begin
			zk := roots[k];

			pz  := PolyEval(coef, zk);
			dpz := PolyDeriv(coef, zk);

			if absC(dpz) < 1e-112 then
			begin
				dpz.re := 1e-110;
				dpz.im := 1e-110;
			end;

			correction_inv := dpz / pz;

			sigma := C(0,0);

			for j := 0 to n-1 do
			begin
				if j = k then Continue;

				diff := zk - roots[j];
				if absC(diff) < 1e-110 then diff.re := 1e-112;

				sigma := sigma + (C(1,0) / diff);
			end;

			correction_inv := correction_inv - sigma;

			if absC(correction_inv) > 1e-112 then
				w := C(1,0) / correction_inv
			else
				w := C(0,0);

			roots[k] := zk - w;

			if absC(w) < 1e-110 then
				Inc(converged);
		end;
	end;

	Writeln('Aberth-Ehrlich iterations: ', iter);

	for k := 0 to n-1 do
		roots_out[k] := roots[k];
end;


var
	coef: array[0..3] of DecFloat;
	roots: array of DecFloatC;
	i: Integer;

begin
	coef[0] := -6;
	coef[1] := 11;
	coef[2] := -6;
	coef[3] := 1;

	SetLength(roots, Length(coef)-1);

	FindRootsAberth(coef, roots);

	Writeln('Roots (Aberth-Ehrlich):');
	for i := 0 to High(roots) do
		writeln(roots[i].re.toString(30),'       ', roots[i].im.toString(30), ' *i');

end.
