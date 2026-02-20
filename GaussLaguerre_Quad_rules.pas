program GaussLaguerreQuadrature;

{$mode objfpc}{$H+}
{$MINFPCONSTPREC 64}
{$MODESWITCH ADVANCEDRECORDS}

uses sysutils, windows, math, unitDecfloat;

{-----------------------------------------------------------
 Symmetric tridiagonal QL algorithm
-----------------------------------------------------------}
procedure tqli(var d: array of decfloat; var e: array of decfloat; 
                var z: array of decfloat; n: Integer);
var
  i, k, l, m, iter: Integer;
  s, r, p, g, f, dd, c, b: decfloat;
begin
  for i := 1 to n-1 do
    e[i-1] := e[i];
  e[n-1] := 0.0;

  for l := 0 to n-1 do
  begin
    iter := 0;
    repeat
      m := l;
      while m < n-1 do
      begin
        dd := Abs(d[m]) + Abs(d[m+1]);
        if Abs(e[m]) + dd = dd then
          Break;
        Inc(m);
      end;

      if m <> l then
      begin
        if iter = 40 then
        begin
          Writeln('Too many iterations');
          Halt;
        end;
        Inc(iter);

        g := (d[l+1] - d[l]) / (2.0 * e[l]);
        r := Sqrt(g*g + 1.0);
        g := d[m] - d[l] + e[l] / (g + Sign(g)*Abs(r));

        s := 1.0;
        c := 1.0;
        p := 0.0;

        for i := m-1 downto l do
        begin
          f := s * e[i];
          b := c * e[i];

          if Abs(f) >= Abs(g) then
          begin
            c := g / f;
            r := Sqrt(c*c + 1.0);
            e[i+1] := f * r;
            s := 1.0 / r;
            c := c * s;
          end
          else
          begin
            s := f / g;
            r := Sqrt(s*s + 1.0);
            e[i+1] := g * r;
            c := 1.0 / r;
            s := s * c;
          end;

          g := d[i+1] - p;
          r := (d[i] - g) * s + 2.0 * c * b;
          p := s * r;
          d[i+1] := g + p;
          g := c * r - b;

          for k := 0 to n-1 do
          begin
            f := z[k*n + i+1];
            z[k*n + i+1] := s * z[k*n + i] + c * f;
            z[k*n + i] := c * z[k*n + i] - s * f;
          end;
        end;

        d[l] := d[l] - p;
        e[l] := g;
        e[m] := 0.0;
      end;
    until m = l;
  end;
end;

{-----------------------------------------------------------
 Generalized Gauss–Laguerre quadrature
-----------------------------------------------------------}
procedure GaussLaguerreGen(n: Integer; alpha: decfloat);
var
  d, e: array of decfloat;
  z: array of decfloat;
  i: Integer;
  gammaFactor: decfloat;
begin
  if alpha <= -1 then
  begin
    Writeln('alpha must be > -1');
    Exit;
  end;

  SetLength(d, n);
  SetLength(e, n);
  SetLength(z, n*n);

  for i := 0 to n-1 do
  begin
    d[i] := 2.0*i + alpha + 1.0;
    if i > 0 then
      e[i] := Sqrt(i * (i + alpha));
  end;

  for i := 0 to n*n-1 do
    z[i] := 0.0;
  
  for i := 0 to n-1 do
    z[i*n + i] := 1.0;

  tqli(d, e, z, n);

  gammaFactor := fpgamma(alpha + 1);

  Writeln;
  Writeln('Generalized Gauss–Laguerre');
  Writeln('n = ', n, ' alpha = ', alpha.toString(2));
  Writeln;
  Writeln('  Node x_i            Weight w_i');
  Writeln('-------------------------------------------');

  for i := 0 to n-1 do
    Writeln(d[i].toString(42), ' ', (gammaFactor * z[i] * z[i]).toString(42));
end;

{-----------------------------------------------------------
 Main program
-----------------------------------------------------------}
var
  n: Integer;
  alpha: decfloat;
  s: ansiString;
begin
  Write('Enter order n: ');
  ReadLn(n);
  Write('Enter alpha (> -1): ');
  ReadLn(s); alpha:=s;

  GaussLaguerreGen(n, alpha);

  WriteLn('Press Enter to continue...');
  ReadLn;
end.
