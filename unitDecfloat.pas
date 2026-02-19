	{$MODE OBJFPC}{$H+}
	{$MINFPCONSTPREC 64}
	{$MODESWITCH ADVANCEDRECORDS}
	///{$CALLING pascal}

	Unit unitDecfloat;

	Interface
	uses sysutils, windows, math, strutils;

	const bias=1073741824;
	const num_dwords=13+1;
	const num_digits=8*(1+num_dwords);
	const Infinity = 1.0 / 0.0;
	const DIVZ_ERR = 1;
	const EXPO_ERR = 2;
	const EXPU_ERR = 3;

const p10_constants : array [0..18] of int64 = (
	1,
	10,
	100,
	1000,
	10000,
	100000,
	1000000,
	10000000,
	100000000,
	1000000000,
	10000000000,
	100000000000,
	1000000000000,
	10000000000000,
	100000000000000,
	1000000000000000,
	10000000000000000,
	100000000000000000,
	1000000000000000000);

	type decfloat = record
		public
		sign : Int32;
		exponent : Int32;
		mantissa : array[0..num_dwords] of Int32;
		function ToString(places:Int32 = num_digits): ansistring;
		function ToStringExp(places:Int32 = num_digits): ansistring;
	end;

	type decfloatc = record
		re, im : decfloat;
		//public
		function ToString(places:Int32 = num_digits): ansistring;
		function ToStringExp(places:Int32 = num_digits): ansistring;
	end;

	var fppi, fpPihalf, fpPIo6, fpPi2, Pi3o2, fpdtor, fpexp1, fpSQR3, fpTMS3, fpln2, fpconst1olog10, sqrt2pi:decfloat;
		gamma_coeff : array [0..70] of decfloat;
		Bernoulli_fac : array [0..75] of decfloat;
		fponec : decfloatc;

	function si2fp( const m : int64; const dwords_in : Int32=num_dwords) : decfloat;
	function fpadd(const x:decfloat; const y:decfloat; const dwords_in:Int32 = num_dwords) : decfloat;
	function ipower(const x:double; const e:Int32) : double;
	function fp2str(const n:decfloat; const places_in : Int32=90) : ansiString;
	procedure gauss_leg_rule(const n:Int32; var x: array of decfloat ; var w: array of decfloat; const dwords_in:Int32 = num_dwords);
	function getExponent(n:decfloat): int32;
	function fpipow(const x:decfloat; const e:int64; const dwords_in:Int32 = num_dwords) : decfloat;

	function FloatToStr(x : decfloat):string; overload;
	procedure Str(x : decfloat; var s:string); overload;

	procedure writelnq(const n:decfloat; const digits : Int32=56);
	procedure writeq(const n:decfloat; const digits : Int32=56);
	operator := (r : Int16) z : decfloat;
	operator < (x : decfloat; y : decfloat) z : boolean;
	operator := (r : ansistring) z : decfloat;
	operator := (r : decfloat) z : ansistring;
	operator := (r : Int32) z : decfloat;
	operator := (r : int64) z : decfloat;
	operator := (r : double) z : decfloat;
	operator + (x : decfloat; y : decfloat) z : decfloat;
	operator - (x : decfloat) z : decfloat;
	operator - (x : decfloat; y : decfloat) z : decfloat;
	operator * (x : decfloat; y : decfloat) z : decfloat;
	operator * (x : decfloat; y : Int32) z : decfloat;
	operator / (x : decfloat; y : decfloat) z : decfloat;
	operator / (x : decfloat; y : Int32) z : decfloat;
	operator <= (x : decfloat; y : decfloat) z : boolean;
	operator > (x : decfloat; y : decfloat) z : boolean;
	operator >= (x : decfloat; y : decfloat) z : boolean;
	operator <> (x : decfloat; y : decfloat) z : boolean;
	operator = (x : decfloat; y : decfloat) z : boolean;
	function fpsign( const x : decfloat): Int32;
	function sqrt(x : decfloat): decfloat;overload;
	function abs(x : decfloat): decfloat;overload;
	function fpSin(x: decfloat): decfloat;
	function fpCos(x: decfloat): decfloat;
	function fpTan(x: decfloat): decfloat;
	function fpArcTan(x_in: decfloat): decfloat;
	function fpArcSin(x: decfloat): decfloat;
	function fpArcCos(x: decfloat): decfloat;

	function sin(x : decfloat): decfloat;overload;
	function cos(x : decfloat): decfloat;overload;
	function tan(x : decfloat): decfloat;overload;
	function arcsin(x : decfloat): decfloat;overload;
	function arccos(x : decfloat): decfloat;overload;
	function arctan(x : decfloat): decfloat;overload;
	function ArcTan2(const Y, X: decfloat): decfloat; overload;
	Function sinh(const x : decfloat) : decfloat; overload;
	Function cosh(const x : decfloat) : decfloat; overload;
	Function tanh(const x : decfloat) : decfloat; overload;
	function hypot(const Y, X: decfloat): decfloat; overload;
	function modulus(const Y, X: decfloat): decfloat; overload; // same as hypot
	function phase(const Y, X: decfloat): decfloat; overload;
	function log(x : decfloat): decfloat;overload;
	function fpexp(const x2:decfloat; const dwords_in:Int32 = num_dwords):decfloat;
	function exp(x : decfloat): decfloat;overload;

//===================================== complex =========================================
	operator := (r : decfloatc) z : ansistring;
//	operator := (r : Int16) z : decfloatc;
//	operator := (r : Int32) z : decfloatc;
//	operator := (r : decfloatc) z : Int32;
//	operator := (r : int64) z : decfloatc;
//	operator := (r : double) z : decfloatc;
	Operator + (x, y : decfloatc) z : decfloatc;
	Operator + (x : decfloatc; y : decfloat) z : decfloatc;
	Operator + (x : decfloat; y : decfloatc) z : decfloatc;
	Operator - (y : decfloatc) z : decfloatc; //negate
	Operator - (x, y : decfloatc) z : decfloatc;
	Operator - (x : decfloatc; y : Double) z : decfloatc;
	Operator - (x : Double; y : decfloatc) z : decfloatc;
	Operator * (x, y : decfloatc) z : decfloatc;
	Operator * (x : decfloatc; y : decfloat) z : decfloatc;
	Operator * (x : decfloat; y : decfloatc) z : decfloatc;
	Operator / (x, y : decfloatc) z : decfloatc;
	Operator / (x : decfloatc; y : Double) z : decfloatc;
	Operator / (x : Double; y : decfloatc) z : decfloatc;
	function Str2c(const s: ansistring) : decfloatc;
	Function csquare(x : decfloatc) : decfloatc;
	Function reciprocal(x :decfloatc) : decfloatc;
	Function cabs(x : decfloatc) : decfloatc;
	Function csqrt(x : decfloatc) : decfloatc;
	Function cexp(x : decfloatc) : decfloatc;
	Function clog(x : decfloatc) : decfloatc;
	Function cpow(Const x, c : decfloatc) : decfloatc;
	Function clogn(Const x, c : decfloatc) : decfloatc;
	Function clog10(Const x : decfloatc) : decfloatc;
	Function csin(const z : decfloatc) : decfloatc;
	Function ccos(const z : decfloatc) : decfloatc;
	Function ctan(const z : decfloatc) : decfloatc;
	Function cArcSin(const x : decfloatc) : decfloatc;
	Function cArcCos(const x : decfloatc) : decfloatc;
	Function cArcTan(const x : decfloatc) : decfloatc;
	Function ccsc(const x : decfloatc) : decfloatc;
	Function csec(const x : decfloatc) : decfloatc;
	Function ccot(const x : decfloatc) : decfloatc;
	Function cArcCsc(const x : decfloatc) : decfloatc;
	Function cArcSec(const x : decfloatc) : decfloatc;
	Function cArcCot(const x : decfloatc) : decfloatc;

	//HYPERBOLICUS I
	Function csinh(const x : decfloatc) : decfloatc;
	Function ccosh(const x : decfloatc) : decfloatc;
	Function ctanh(const x : decfloatc) : decfloatc;
	Function cArcSinh(const x : decfloatc) : decfloatc;
	Function cArcCosh(const x : decfloatc) : decfloatc;
	Function cArcTanh(const x : decfloatc) : decfloatc;

	//HYPERBOLICUS II
	Function ccsch(const x : decfloatc) : decfloatc;
	Function csech(const x : decfloatc) : decfloatc;
	Function ccoth(const x : decfloatc) : decfloatc;
	Function cArcCsch(const x : decfloatc) : decfloatc;
	Function cArcSech(const x : decfloatc) : decfloatc;
	Function cArcCoth(const x : decfloatc) : decfloatc;
	Function fpsignc(const x : decfloatc) : Int32;

	Operator ** (const x : decfloatc; const y : decfloatc) : decfloatc;
	Operator ** (const x : Double; const y : decfloatc) : decfloatc;
	Operator ** (const x : Integer; const y : decfloatc) : decfloatc;
	function toDecfloatc(const x : decfloat; const y : decfloat) : decfloatc; overload;
	function toDecfloatc(const x : double; const y : decfloat) : decfloatc; overload;
	function toDecfloatc(const x : decfloat; const y : double) : decfloatc; overload;
	function toDecfloatc(const x : double; const y : double) : decfloatc; overload;
	function toDecfloatc(const x : Int32; const y : decfloat) : decfloatc; overload;
	function toDecfloatc(const x : decfloat; const y : Int32) : decfloatc; overload;
	function toDecfloatc(const x : Int32; const y : Int32) : decfloatc; overload;
	function cArctan2(y, x: decfloatc): decfloatc; overload;

	function fpgamma(x : decfloat) : decfloat;
	function fpgammac(x : decfloatc) : decfloatc;
	Function fpzeta(n : decfloat) : decfloat;
	Function fpzetac(n : decfloatc) : decfloatc;
	function log_agm(const x:decfloat; const digits:Int32=num_digits) : decfloat;
	operator := (r : decfloat) z : Int32;
	function frac(x : decfloat): decfloat;overload;
	function trunc(x : decfloat): decfloat;overload;

	implementation

	procedure copydf(var resul:decfloat; const source:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;
		resul.sign := source.sign;
		resul.exponent := source.exponent;
		for i := 0 to dwords do resul.mantissa[i] := source.mantissa[i];
	end;

	function stringn(const n : Int32; const c:char ): ansistring;
	var
		i:Int32;
		s:ansistring;
	begin
		s:='';
		for i:=1 to n do
		begin
			s:=s+c;
		end;
		result:=s;
	end;

	function log10_32(const value : longword) : Int32;
	begin
		result:=-1;
		if (value >= 1000000000) then exit(9);
		if (value >= 100000000) then exit(8);
		if (value >= 10000000) then exit(7);
		if (value >= 1000000) then exit(6);
		if (value >= 100000) then exit(5);
		if (value >= 10000) then exit(4);
		if (value >= 1000) then exit(3);
		if (value >= 100) then exit(2);
		if (value >= 10) then exit(1);
		if (value >= 1) then exit(0);
	end;

	function log10_64(const value : qword) : Int32;
	begin
		result:=-1;
		if (value >= 10000000000000000000) then exit(19);
		if (value >= 1000000000000000000) then exit(18);
		if (value >= 100000000000000000) then exit(17);
		if (value >= 10000000000000000) then exit(16);
		if (value >= 1000000000000000) then exit(15);
		if (value >= 100000000000000) then exit(14);
		if (value >= 10000000000000) then exit(13);
		if (value >= 1000000000000) then exit(12);
		if (value >= 100000000000) then exit(11);
		if (value >= 10000000000) then exit(10);
		if (value >= 1000000000) then exit(9);
		if (value >= 100000000) then exit(8);
		if (value >= 10000000) then exit(7);
		if (value >= 1000000) then exit(6);
		if (value >= 100000) then exit(5);
		if (value >= 10000) then exit(4);
		if (value >= 1000) then exit(3);
		if (value >= 100) then exit(2);
		if (value >= 10) then exit(1);
		if (value >= 1) then exit(0);
	end;

	function fp2str_exp (const n : decfloat; const places_in : Int32=num_digits) : ansiString;
	var
		places, i, ex : Int32;
		v, f, ts :ansistring;
	begin
		places:=places_in;
		If (places>num_digits) Then places:=num_digits;

		If n.exponent > 0 Then
			ex := (n.exponent and $7fffffff) - bias - 1
		Else
			ex := 0;

		If n.sign<>0 Then
			v := '-'
		Else
			v := ' ';

		Str(n.mantissa[0], ts);
		ts:=trim(ts);

		while Length(ts) < 8 do
		begin
			ts := ts + '0';
		End;
		v := v + ts[1]+'.'+copy(ts, 2, length(ts)-1);

		For i := 1 To num_dwords do
		begin
			Str(n.mantissa[i], ts);
			ts:=trim(ts);
			while Length(ts) < 8 do
			begin
				ts := '0' + ts;
			End;
			v := v + ts;
		end;

		v:=copy(v, 1, places + 2);
		v:=v.trimright(['0']);
		Str(Abs(ex), f);
		f:=trim(f);

		If ex < 0 Then
			v := v + 'e-'
		Else
			v := v + 'e+';

		v := v + f;

		result:=v;
	end;

	function fpfix(const num:decfloat; const dwords_in:Int32 = num_dwords) : decfloat;
	var
		dwords, ex, ex2, j, k:Int32;
		ip:decfloat;

	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for j:=0 to num_dwords do
			ip.mantissa[j]:=0;

		ex := (num.exponent and $7FFFFFFF) - bias;
		if (ex < 1) then
		begin
			result:=si2fp(0);
			exit(result);
		end;
		if (ex >= (8 * dwords)) then exit(num);
		ex2 := ex div 8;
		k := ex2;
		j := ex mod 8;
		while (ex2 > 0) do
		begin
			ex2 -= 1;
			ip.mantissa[ex2] := num.mantissa[ex2];
		end;
		if (j = 1) then
			ip.mantissa[k] := 10000000 * (num.mantissa[k] div 10000000)
		else if (j = 2) then
			ip.mantissa[k] := 1000000 * (num.mantissa[k] div 1000000)
		else if (j = 3) then
			ip.mantissa[k] := 100000 * (num.mantissa[k] div 100000)
		else if (j = 4) then
			ip.mantissa[k] := 10000 * (num.mantissa[k] div 10000)
		else if (j = 5) then
			ip.mantissa[k] := 1000 * (num.mantissa[k] div 1000)
		else if (j = 6) then
			ip.mantissa[k] := 100 * (num.mantissa[k] div 100)
		else if (j = 7) then
			ip.mantissa[k] := 10 * (num.mantissa[k] div 10)
		else if (j = 8) then
			ip.mantissa[k] := num.mantissa[k];

		ip.exponent := ex + bias;
		ip.sign := num.sign;

		result := ip;
	end;

	function fpfix_exp(const num:decfloat; const places_in:Int32 = num_digits) : ansistring;
	var
		places, ex, ex2, j, k, d:Int32;
		ip, tmp:decfloat;

	begin
		places := places_in;
		d:=0;
		if (places < (num_digits-2)) then places:=places+1;

		for j:=0 to num_dwords do
			ip.mantissa[j]:=0;

		ex := places;
		if ex<1 then
			exit ('0');

		if ex>=(num_digits) then
			exit (fp2str_exp (num, places));

		ex2:=ex div 8;
		k:=ex2;
		j:=ex mod 8;
		while ex2>0 do
		begin
			ex2:=ex2-1;
			ip.mantissa[ex2]:=num.mantissa[ex2];
		end;

		if j=0 then
			d:=num.mantissa[k] div 10000000
		else if j=1 then
		begin
			ip.mantissa[k]:=10000000*(num.mantissa[k] div 10000000);
			d:=(num.mantissa[k] div 1000000) mod 10;
		end
		else if j=2 then
		begin
			ip.mantissa[k]:=1000000*(num.mantissa[k] div 1000000);
			d:=(num.mantissa[k] div 100000) mod 10;
		end
		else if j=3 then
		begin
			ip.mantissa[k]:=100000*(num.mantissa[k] div 100000);
			d:=(num.mantissa[k] div 10000) mod 10;
		end
		else if j=4 then
		begin
			ip.mantissa[k]:=10000*(num.mantissa[k] div 10000);
			d:=(num.mantissa[k] div 1000) mod 10;
		end
		else if j=5 then
		begin
			ip.mantissa[k]:=1000*(num.mantissa[k] div 1000);
			d:=(num.mantissa[k] div 100) mod 10;
		end
		else if j=6 then
		begin
			ip.mantissa[k]:=100*(num.mantissa[k] div 100);
			d:=(num.mantissa[k] div 10) mod 10;
		end
		else if j=7 then
		begin
			ip.mantissa[k]:=10*(num.mantissa[k] div 10);
			d:=(num.mantissa[k]) mod 10;
		end;
		ip.exponent:=ex+bias;
		tmp:=si2fp(1);
		if d>=5 then
			ip:=fpadd(ip, tmp);

		if ip.exponent>ex+bias then
			ip.exponent:=num.exponent+1
		else
			ip.exponent:=num.exponent;

		ip.sign:=num.sign;
		Result:= fp2str_exp (ip, places-1);
	end;

	function fpfix_is_odd(const num:decfloat) : Int32;
	var
		ex, ex2, j, k: Int32;
	begin
		ex := (num.exponent and $7fffffff) - bias;
		if (ex < 1) then exit(0);
		if (ex >= num_digits) then exit(-1);
		ex2 := ex div 8;
		k := ex2;
		j := ex mod 8;

		if (j = 1) then exit((num.mantissa[k] div 10000000) and 1);
		if (j = 2) then exit((num.mantissa[k] div 1000000) and 1);
		if (j = 3) then exit((num.mantissa[k] div 100000) and 1);
		if (j = 4) then exit((num.mantissa[k] div 10000) and 1);
		if (j = 5) then exit((num.mantissa[k] div 1000) and 1);
		if (j = 6) then exit((num.mantissa[k] div 100) and 1);
		if (j = 7) then exit((num.mantissa[k] div 10) and 1);
		if (j = 8) then exit(num.mantissa[k] and 1);

		result := 0;
	end;

	function fp2dbl(const n:decfloat) : double;
	var
		dbl, tp : double;
		ex, sign : Int32;
	begin
		if (n.sign <> 0) then
			sign := -1
		else
			sign := 1;

		if (n.exponent > 0) then
			ex := (n.exponent and $7fffffff) - bias - 1
		else
			ex := 0;

		if (ex > 308) then exit(infinity);
		if (ex < (-308)) then exit(0);
		dbl := double(n.mantissa[0]) + (double(n.mantissa[1]) / 100000000.0);
		dbl := dbl / 10000000.0;
		tp := ipower(10.0, ex);
		result := dbl * tp * sign;
	end;

	function str2fp(const value_in : ansistring) : decfloat;
	var
		j, s, d, e, ep, ulng  : Int32;
		ex, es, i, f, fp, fln, cv, cc  : Int32;
		f1, f2, f3, ts, value : ansistring;
		c : string;
		n : decfloat;
	begin
		value:=UpperCase(value_in.trimleft(['0']));
		if length(value)=0 then value:='0';
		fln:=length(value);
		n.sign:=0;
		n.exponent:=0;
		j := 1;
		s := 1;
		d := 0;
		e := 0;
		ep := 0;
		ex := 0;
		es := 1;
		i := 1;
		f := 0;
		fp := 0;
		f1 := '';
		f2 := '';
		f3 := '';
		c:= '';
		cc:=0;

		while (j<=fln) do
		begin
			c:=value[j];
			if (ep=1) then
			begin
				if (c=' ') then
				begin
					inc(j);
					continue;
				end;
				if (c='-') then
				begin
					es:=-es;
					c:='';
					inc(j);
					continue;
				end;
				if (c='+') then
				begin
					inc(j);
					continue;
				end;
				if ((c='0') and (length(f3)=0)) then
				begin
					inc(j);
					continue;
				end;
				if ((c>'/') and (c<':')) then
				begin
					f3:=f3+c;
					val(c, cv, cc);
					ex:=10*ex+cv;
					inc(j);
					continue;
				end;
			end;
			if (c=' ') then
			begin
				inc(j);
				continue;
			end;
			if (c='-') then
			begin
				s:=-s;
				inc(j);
				continue;
			end;
			if (c='+') then
			begin
				inc(j);
				continue;
			end;
			if (c='.') then
			begin
				if (d=1) then
				begin
					inc(j);
					continue;
				end;
				d:=1;
			end;
			if ((c>'/') and (c<':')) then
			begin
				If ((c = '0') And (i = 0)) Then
				begin
					if (d=0) then
					begin
						inc(j);
						continue;
					end;
					if ((d=1) and (f=0))then
					begin
						dec(e);
						inc(j);
						continue;
					end;
				end;
				if (d=0) then
				begin
					f1:=f1+c;
					inc(i);
				end;
				if (d<>0) then
				begin
					if (c>'0') then fp:=1;
					f2:=f2+c;
					inc(f);
				end;
			end;
			if (c='E') then
			begin
				ep:=1;
			end;
			inc(j);
		end;
		if (fp=0) then
		begin
			f:=0;
			f2:='';
		end;
		if (s=1) then s:=0;
		n.sign:=cc;
		n.sign:=s;
		ex := es * ex - 2 + i + e;

		f1 := f1 + f2;
		fln:=length(f1);
		if (fln>num_digits) then
		begin
			f1:=copy(f1, 1, num_digits);
		end;
		while (length(f1)<num_digits) do
		begin
			f1:=f1+'0';
		end;
		j:=1;
		for i:=0 to num_dwords do
		begin
			ts:=copy(f1, j, 8);
			val(ts, ulng, cc);
			n.mantissa[i]:=ulng;
			if (ulng<>0) then fp:=1;
			j:=j+8;
		end;

		if (fp>0) then n.exponent:=ex+bias+1;
		result:=n;
	end;

	function si2fp( const m : int64; const dwords_in : Int32=num_dwords) : decfloat;
	var
		dwords, ex, i : Int32;
		mm : int64;
	begin
		dwords:=dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;
		mm := abs(m);

		for i := 0 to dwords do
		begin
			result.mantissa[i] := 0;
		end;

		if (m = 0) then
		begin
			result.exponent := 0;
			result.sign := 0;
			exit(result);
		end;
		ex := (log10_64(mm));
		result.exponent := bias + ex + 1;

		if (ex > 15) then
		begin
			ex := ex-7;
			result.mantissa[0] := (mm div int64(p10_constants[ex]));
			mm := mm mod int64(p10_constants[ex]);
			ex := ex-8;
			result.mantissa[1] := (mm div int64(p10_constants[ex]));
			mm := mm mod int64(p10_constants[ex]);
			result.mantissa[2] := (mm * int64(p10_constants[8 - ex]));
		end
		else if (ex >= 8) then
		begin
			ex := ex-7;
			result.mantissa[0] := (mm div int64(p10_constants[ex]));
			mm := mm mod int64(p10_constants[ex]);
			result.mantissa[1] := (mm * int64(p10_constants[8 - ex]));
		end
		else if (ex < 8) then
		begin
			result.mantissa[0] := (mm * int64(p10_constants[7 - ex]));
		end;
		result.sign := 0;
		if (m < 0) then result.sign := -1;
	end;

	procedure rshift_1(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := dwords downto 1 do
		begin
			v1 := mantissa.mantissa[i] div 10;
			v2 := mantissa.mantissa[i - 1] mod 10;
			v2 := v2 * 10000000 + v1;
			mantissa.mantissa[i] := v2;
		end;
		mantissa.mantissa[0] := mantissa.mantissa[0] div 10;
	end;

	procedure lshift_1(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := 0 to (dwords - 1) do
		begin
			v1 := mantissa.mantissa[i] mod 10000000;
			v2 := mantissa.mantissa[i + 1] div 10000000;
			mantissa.mantissa[i] := v1 * 10 + v2;
			mantissa.mantissa[i + 1] := mantissa.mantissa[i + 1] mod 10000000;
		end;
		mantissa.mantissa[dwords] := 10 * (mantissa.mantissa[dwords] mod 10000000);
	end;

	procedure rshift_2(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := dwords downto 1 do
		begin
			v1 := mantissa.mantissa[i] div 100;
			v2 := mantissa.mantissa[i - 1] mod 100;
			v2 := v2 * 1000000 + v1;
			mantissa.mantissa[i] := v2;
		end;
		mantissa.mantissa[0] := mantissa.mantissa[0] div 100;
	end;

	procedure lshift_2(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := 0 to (dwords - 1) do
		begin
			v1 := mantissa.mantissa[i] mod 1000000;
			v2 := mantissa.mantissa[i + 1] div 1000000;
			mantissa.mantissa[i] := v1 * 100 + v2;
			mantissa.mantissa[i + 1] := mantissa.mantissa[i + 1] mod 1000000;
		end;
		mantissa.mantissa[dwords] := 100 * (mantissa.mantissa[dwords] mod 1000000);
	end;

	procedure rshift_3(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := dwords downto 1 do
		begin
			v1 := mantissa.mantissa[i] div 1000;
			v2 := mantissa.mantissa[i - 1] mod 1000;
			v2 := v2 * 100000 + v1;
			mantissa.mantissa[i] := v2;
		end;
		mantissa.mantissa[0] := mantissa.mantissa[0] div 1000;
	end;

	procedure lshift_3(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := 0 to (dwords - 1) do
		begin
			v1 := mantissa.mantissa[i] mod 100000;
			v2 := mantissa.mantissa[i + 1] div 100000;
			mantissa.mantissa[i] := v1 * 1000 + v2;
			mantissa.mantissa[i + 1] := mantissa.mantissa[i + 1] mod 100000;
		end;
		mantissa.mantissa[dwords] := 1000 * (mantissa.mantissa[dwords] mod 100000);
	end;

	procedure rshift_4(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := dwords downto 1 do
		begin
			v1 := mantissa.mantissa[i] div 10000;
			v2 := mantissa.mantissa[i - 1] mod 10000;
			v2 := v2 * 10000 + v1;
			mantissa.mantissa[i] := v2;
		end;
		mantissa.mantissa[0] := mantissa.mantissa[0] div 10000;
	end;

	procedure lshift_4(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := 0 to (dwords - 1) do
		begin
			v1 := mantissa.mantissa[i] mod 10000;
			v2 := mantissa.mantissa[i + 1] div 10000;
			mantissa.mantissa[i] := v1 * 10000 + v2;
			mantissa.mantissa[i + 1] := mantissa.mantissa[i + 1] mod 10000;
		end;
		mantissa.mantissa[dwords] := 10000 * (mantissa.mantissa[dwords] mod 10000);
	end;

	procedure rshift_5(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := dwords downto 1 do
		begin
			v1 := mantissa.mantissa[i] div 100000;
			v2 := mantissa.mantissa[i - 1] mod 100000;
			v2 := v2 * 1000 + v1;
			mantissa.mantissa[i] := v2;
		end;
		mantissa.mantissa[0] := mantissa.mantissa[0] div 100000;
	end;

	procedure lshift_5(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := 0 to (dwords - 1) do
		begin
			v1 := mantissa.mantissa[i] mod 1000;
			v2 := mantissa.mantissa[i + 1] div 1000;
			mantissa.mantissa[i] := v1 * 100000 + v2;
			mantissa.mantissa[i + 1] := mantissa.mantissa[i + 1] mod 1000;
		end;
		mantissa.mantissa[dwords] := 100000 * (mantissa.mantissa[dwords] mod 1000);
	end;

	procedure rshift_6(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := dwords downto 1 do
		begin
			v1 := mantissa.mantissa[i] div 1000000;
			v2 := mantissa.mantissa[i - 1] mod 1000000;
			v2 := v2 * 100 + v1;
			mantissa.mantissa[i] := v2;
		end;
		mantissa.mantissa[0] := mantissa.mantissa[0] div 1000000;
	end;

	procedure lshift_6(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := 0 to (dwords - 1) do
		begin
			v1 := mantissa.mantissa[i] mod 100;
			v2 := mantissa.mantissa[i + 1] div 100;
			mantissa.mantissa[i] := v1 * 1000000 + v2;
			mantissa.mantissa[i + 1] := mantissa.mantissa[i + 1] mod 100;
		end;
		mantissa.mantissa[dwords] := 1000000 * (mantissa.mantissa[dwords] mod 100);
	end;

	procedure rshift_7(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := dwords downto 1 do
		begin
			v1 := mantissa.mantissa[i] div 10000000;
			v2 := mantissa.mantissa[i - 1] mod 10000000;
			v2 := v2 * 10 + v1;
			mantissa.mantissa[i] := v2;
		end;
		mantissa.mantissa[0] := mantissa.mantissa[0] div 10000000;
	end;

	procedure lshift_7(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v1, v2:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := 0 to (dwords - 1) do
		begin
			v1 := mantissa.mantissa[i] mod 10;
			v2 := mantissa.mantissa[i + 1] div 10;
			mantissa.mantissa[i] := v1 * 10000000 + v2;
			mantissa.mantissa[i + 1] := mantissa.mantissa[i + 1] mod 10;
		end;
		mantissa.mantissa[dwords] := 10000000 * (mantissa.mantissa[dwords] mod 10);
	end;

	procedure rshift_8(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := dwords downto 1 do
		begin
			mantissa.mantissa[i] := mantissa.mantissa[i - 1];
		end;
		mantissa.mantissa[0] := 0;
	end;

	procedure lshift_8(var mantissa:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		for i := 0 to (dwords - 1) do
		begin
			mantissa.mantissa[i] := mantissa.mantissa[i + 1];
		end;
		mantissa.mantissa[dwords] := 0;
	end;

	function fpcmp(const x:decfloat; const y:decfloat; const dwords_in:Int32=num_dwords) :Int32;
	var
		dwords, i:Int32;
		c:int64;
	begin
		dwords := dwords_in;
		c:=0;
		if (dwords > num_dwords) then dwords := num_dwords;
		if (x.sign < y.sign) then exit(-1);
		if (x.sign > y.sign) then exit(1);
		if (x.exponent < y.exponent) then
		begin
			if (x.sign = 0) then
				exit(-1)
			else
				exit(1);
		end;
		if (x.exponent > y.exponent) then
		begin
			if (x.sign = 0) then
				exit(1)
			else
				exit(-1);
		end;

		for i := 0 to dwords do
		begin
			c := int64(x.mantissa[i]) - int64(y.mantissa[i]);
			if (c <> 0) then break;
		end;
		if (c = 0) then exit(0);
		if (c < 0) then
		begin
			if (x.sign = 0) then
				exit(-1)
			else
				exit(1);
		end;
		if (c > 0) then
		begin
			if (x.sign = 0) then
				exit(1)
			else
				exit(-1);
		end;
		exit(0);
	end;

	function fpcmp_abs(const x:decfloat; const y:decfloat; const dwords_in:Int32 = num_dwords) :Int32;
	var
		dwords, i:Int32;
		c:int64;
	begin
		dwords := dwords_in;
		c:=0;
		if (dwords > num_dwords) then dwords := num_dwords;

		if (x.exponent < y.exponent) then exit(-1);
		if (x.exponent > y.exponent) then exit(1);

		for i := 0 to dwords do
		begin
			c := int64(x.mantissa[i]) - int64(y.mantissa[i]);
			if (c <> 0) then break;
		end;
		if (c = 0) then exit(0);
		if (c < 0) then exit(-1);
		if (c > 0) then exit(1);
		exit(0);
	end;

	function norm_fac1(var fac1:decfloat; const dwords_in:Int32 = num_dwords) : Int32;
	var
		dwords, i, f:Int32;
	begin
		dwords := dwords_in;
		f:=0;
		if (dwords > num_dwords) then dwords := num_dwords;
		// normalize the number in fac1
		// all routines exit through this one.

		//see if the mantissa is all zeros.
		//if so, set the exponent and sign equal to 0.

		for i := 0 to dwords do
		begin
			if (fac1.mantissa[i] > 0) then f := 1;
		end;
		if (f = 0) then
		begin
			fac1.exponent := 0;
			fac1.sign := 0;
			exit(0);
			//if the highmost digit in fac1_man is nonzero,
			//shift the mantissa right 1 digit and
			//increment the exponent
		end
		else if (fac1.mantissa[0] > 99999999) then
		begin
			rshift_1(fac1, dwords);
			inc(fac1.exponent);
		end
		else
		begin
			//now shift fac1_man 1 to the left until a
			//nonzero digit appears in the }-to-highest
			//digit of fac1_man.  decrement exponent for
			//each shift.
			while (fac1.mantissa[0] = 0) do
			begin
				lshift_8(fac1, dwords);
				dec(fac1.exponent, 8);
				if (fac1.exponent = 0) then
					exit(EXPU_ERR);
			end;
			if (fac1.mantissa[0] < 10) then
			begin
				lshift_7(fac1, dwords);
				dec(fac1.exponent, 7);
			end
			else if (fac1.mantissa[0] < 100) then
			begin
				lshift_6(fac1, dwords);
				dec(fac1.exponent, 6);
			end
			else if (fac1.mantissa[0] < 1000) then
			begin
				lshift_5(fac1, dwords);
				dec(fac1.exponent, 5);
			end
			else if (fac1.mantissa[0] < 10000) then
			begin
				lshift_4(fac1, dwords);
				dec(fac1.exponent, 4);
			end
			else if (fac1.mantissa[0] < 100000) then
			begin
				lshift_3(fac1, dwords);
				dec(fac1.exponent, 3);
			end
			else if (fac1.mantissa[0] < 1000000) then
			begin
				lshift_2(fac1, dwords);
				dec(fac1.exponent, 2);
			end
			else if (fac1.mantissa[0] < 10000000) then
			begin
				lshift_1(fac1, dwords);
				dec(fac1.exponent);
			end;
		end;
		//check for overflow/underflow
		if (fac1.exponent < 0) then
			exit(EXPO_ERR);

		exit(0);
	end;

	procedure fpadd_aux(var fac1:decfloat; var fac2:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v, c:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		c := 0;
		for i := dwords downto 1 do
		begin
			v := fac2.mantissa[i] + fac1.mantissa[i] + c;
			if (v > 99999999) then
			begin
				v := v - 100000000;
				c := 1;
			end
			else
				c := 0;

			fac1.mantissa[i] := v;
		end;
		v := fac1.mantissa[0] + fac2.mantissa[0] + c;
		fac1.mantissa[0] := v;
		i := norm_fac1(fac1, dwords);
	end;

	procedure fpsub_aux(var fac1:decfloat; var fac2:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords, i, v, c:Int32;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		c := 0;
		for i := dwords downto 1 do
		begin
			v := fac1.mantissa[i] - fac2.mantissa[i] - c;
			if (v < 0) then
			begin
				v := v + 100000000;
				c := 1;
			end
			else
				c := 0;

			fac1.mantissa[i] := v;
		end;
		v := fac1.mantissa[0] - fac2.mantissa[0] - c;
		fac1.mantissa[0] := v;
		i := norm_fac1(fac1, dwords);
	end;

	function fpadd(const x:decfloat; const y:decfloat; const dwords_in:Int32 = num_dwords) : decfloat;
	var
		dwords, i, t, c:Int32;
		fac1, fac2:decfloat;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		c := fpcmp_abs(x, y, dwords);

		if (c < 0) then
		begin
			fac1 := y;
			fac2 := x;
		end
		else
		begin
			fac1 := x;
			fac2 := y;
		end;
		t := fac1.exponent - fac2.exponent;

		t := ((fac1.exponent and $7fffffff) - bias - 1) - ((fac2.exponent and $7fffffff) - bias - 1);

		if (t < (num_digits+8)) then
		begin
			//the difference between the two
			//exponents indicate how many times
			//we have to multiply the mantissa
			//of fac2 by 10 (i.e., shift it right 1 place).
			//if we have to shift more times than
			//we have dwords, the result is already in fac1.
			t := fac1.exponent - fac2.exponent;
			if ((t > 0) and (t < (num_digits+8))) then //shift
			begin
				i := t div 8;
				while (i > 0) do
				begin
					rshift_8(fac2, dwords);
					dec(t, 8);
					dec(i);
				end;
				if (t = 7) then rshift_7(fac2, dwords)
				else if (t = 6) then rshift_6(fac2, dwords)
				else if (t = 5) then rshift_5(fac2, dwords)
				else if (t = 4) then rshift_4(fac2, dwords)
				else if (t = 3) then rshift_3(fac2, dwords)
				else if (t = 2) then rshift_2(fac2, dwords)
				else if (t = 1) then rshift_1(fac2, dwords);
			end;
			//see if the signs of the two numbers
			//are the same.  if so, add; if not, subtract.
			if (fac1.sign = fac2.sign) then //add
			begin
				fpadd_aux(fac1, fac2, dwords);
			end
			else
			begin
				fpsub_aux(fac1, fac2, dwords)
			end;
		end;
		result := fac1;
	end;

	function fpsub(const x:decfloat; const y:decfloat; const dwords_in:Int32 = num_dwords) : decfloat;
	var
		dwords:Int32;
		fac1, fac2:decfloat;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		fac1 := x;
		fac2 := y;
		fac2.sign := (fac2.sign xor (-1));
		result := fpadd(fac1, fac2, dwords);
	end;

	procedure LSHIFT_a (var mantissa:array of UInt32; const n:Int32; const k:Int32);
		var
		v1, v2, c1, c2 :UInt32;
		i :Int32;
		begin
			c1 := p10_constants[k];
			c2 := p10_constants[8 - k];
			For i := 0 To n do
			begin
				v1 := mantissa[i] mod c2;
				v2 := mantissa[i + 1] div c2;
				mantissa[i] := (v1 * c1 + v2);
				mantissa[i + 1] := mantissa[i + 1] mod c2;
			end;
			mantissa[n+1] := c1 * (mantissa[n+1] mod c2);
	End;

	function fpmul(const x:decfloat; const y:decfloat; const dwords_in:Int32 = num_dwords) : decfloat;
	var
		dwords, i, j, ex, den, num:Int32;
		digit, carry, prod:int64;
		fac1, fac2:decfloat;
		fac3 : array of longword;
		//fac3 : array[0..2 * num_dwords + 2] of longword;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;
		setLength(fac3, 2 * num_dwords + 3);

		//FillChar(result.mantissa[0], SizeOf(result.mantissa[0]) * Length(result.mantissa), 0);
		result:=si2fp(0);
		ex := 0;
		carry := 0;
		//uint32_t fac3[dwords*2+2];

		fac1 := x;
		fac2 := y;

		//check exponents.  if either is zero,
		//the result is zero
		if ((fac1.exponent = 0) or (fac2.exponent = 0)) then //result is zero...clear fac1.
		begin
			fac1.sign := 0;
			fac1.exponent := 0;
			for i := 0 to dwords do fac1.mantissa[i] := 0;
			//norm_fac1(fac1)
			exit(fac1);
		end
		else
		begin
			if (ex < 0) then
			begin
				//er:=EXPO_ERR;
				exit(fac1);
			end;

			//clear fac3 mantissa
			for i := 0 to (dwords * 2 + 2) do fac3[i] := 0;

			den := dwords;
			while (fac2.mantissa[den] = 0) do dec(den);
			num := dwords;
			while (fac1.mantissa[num] = 0) do dec(num);

			if (num < den) then
			begin
				fac1:=y;
				fac2:=x;
				i := den; den := num;
				num := i;
			end;

			for j := den downto 0 do
			begin
				carry := 0;
				digit := int64(fac2.mantissa[j]);
				for i := num downto 0 do
				begin
					prod := int64(fac3[i + j + 1]) + digit * int64(fac1.mantissa[i]) + carry;
					carry := prod div 100000000;
					fac3[i + j + 1] := longword(prod mod 100000000);
				end;

				fac3[j] := longword(carry);
			end;

			for i := 0 to dwords do fac1.mantissa[i] := fac3[i];
		end;

		j := 0;
		If carry=0 Then
		begin
			j := 8;
			For i := 0 To dwords+2 do fac3[i] := fac3[i+1];
		End;

		if fac3[0] < 10 then j := 7
		else if fac3[0] < 100 then j := 6
		else if fac3[0] < 1000 then j := 5
		else if fac3[0] < 10000 then j := 4
		else if fac3[0] < 100000 then j := 3
		else if fac3[0] < 1000000 then j := 2
		else if fac3[0] < 10000000 then j := 1;

		If j>0 Then LSHIFT_a(fac3, dwords+j, j);

		For i := 0 To dwords do result.mantissa[i] := fac3[i];

		(* now determine exponent of result. *)
		(* as you do...watch for overflow. *)
		ex := fac2.exponent-bias+fac1.exponent;
		result.exponent := ex-j;
		(* determine the sign of the product *)
		result.sign := fac1.sign Xor fac2.sign;

		norm_fac1(result, dwords);
	end;

	function fpmul_si(const x:decfloat; const y:Int32; const dwords_in:Int32 = num_dwords) : decfloat;
	var
		dwords, i:Int32;
		carry, digit, prod, value:int64;
		fac1, fac2:decfloat;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		fac1 := x;

		if (y < 0) then
			digit := -y
		else
			digit := y;

		if (digit > 99999999) then
		begin
			fac2:=si2fp(y, dwords);
			result:=fpmul(fac1, fac2, dwords);
			exit(result);
		end;
		//check exponents.  if either is zero,
		//the result is zero
		if ((fac1.exponent = 0) or (y = 0)) then //result is zero...clear fac1.
		begin
			result:=si2fp(0, dwords);
			exit(result);
		end;

		if (digit = 1) and (y < 0)then
		begin
			fac1.sign := fac1.sign xor (-1);
			exit(fac1);
		end;

		carry := 0;

		for i := dwords downto 0 do
		begin
			prod := digit * fac1.mantissa[i] + carry;
			value := (prod mod 100000000);
			fac1.mantissa[i] := Int32(value);
			carry := prod div 100000000;
		end;

		if (carry < 10) then
		begin
			rshift_1(fac1, dwords);
			fac1.exponent += 1;
			fac1.mantissa[0] += Int32(carry) * 10000000;
		end
		else if (carry < 100) then
		begin
			rshift_2(fac1, dwords);
			fac1.exponent += 2;
			fac1.mantissa[0] += Int32(carry) * 1000000;
		end
		else if (carry < 1000) then
		begin
			rshift_3(fac1, dwords);
			fac1.exponent += 3;
			fac1.mantissa[0] += Int32(carry) * 100000;
		end
		else if (carry < 10000) then
		begin
			rshift_4(fac1, dwords);
			fac1.exponent += 4;
			fac1.mantissa[0] += Int32(carry) * 10000;
		end
		else if (carry < 100000) then
		begin
			rshift_5(fac1, dwords);
			fac1.exponent += 5;
			fac1.mantissa[0] += Int32(carry) * 1000;
		end
		else if (carry < 1000000) then
		begin
			rshift_6(fac1, dwords);
			fac1.exponent += 6;
			fac1.mantissa[0] += Int32(carry) * 100;
		end
		else if (carry < 10000000) then
		begin
			rshift_7(fac1, dwords);
			fac1.exponent += 7;
			fac1.mantissa[0] += Int32(carry) * 10;
		end
		else if (carry < 100000000) then
		begin
			rshift_8(fac1, dwords);
			fac1.exponent += 8;
			fac1.mantissa[0] += Int32(carry);
		end;

		i := norm_fac1(fac1, dwords);

		if (y < 0) then
			fac1.sign := fac1.sign xor (-1);
		result := fac1;
	end;

	procedure LSHIFT_da(var mantissa:array of double; const n:Int32; const k:Int32);
	var
		v1, v2, c1, c2:longword;
		i:Int32;
	begin
		c1 := longword(p10_constants[k]);
		c2 := longword(p10_constants[4 - k]);
		for i := 2 to (n - 1) do
		begin
			v1 := round(mantissa[i]) mod c2;
			v2 := round(mantissa[i + 1]) div c2;
			mantissa[i] := v1 * c1 + v2;
			mantissa[i + 1] := round(mantissa[i + 1]) mod c2;
		end;
		mantissa[n] := c1 * (round(mantissa[n]) mod c2);
	end;

	function realw(var w:array of double; const j:Int32; const ubw:Int32) : double;
	var
		wx:double;
	begin
		wx := ((w[j - 1] * 10000 + w[j]) * 10000 + w[j + 1]) * 10000;
		if (ubw >= (j + 2)) then
			wx := wx + w[j + 2];

		result:=wx;
	end;

	procedure subtract(var w:array of double; const q:Int32; var d:array of double; const ka:Int32; const kb:Int32);
	var
		j:Int32;
	begin
		for j := ka to kb do
			w[j] := w[j] - q * d[j - ka + 2];
	end;

	procedure normalize(var w:array of double; const ka:Int32; const q:Int32);
	begin
		w[ka] := w[ka] + w[ka - 1] * 10000;
		w[ka - 1] := q;
	end;

	procedure finalnorm(var w:array of double; const kb:Int32);
	const b = 10000;
	var
		carry, j:Int32;
	begin
		for j := kb downto 3 do
		begin
			if w[j] < 0 then
				carry := round((-w[j] - 1) / b) + 1
			else
			begin
				if w[j] >= b then
					carry := -round(w[j] / b)
				else
					carry := 0;
			end;
			w[j] := w[j] + carry * b;
			w[j - 1] := w[j - 1] - carry;
		end;
	end;

	function fpdiv(const x:decfloat; const y:decfloat; const dwords_in:Int32 = num_dwords) : decfloat;

	const b = 10000;
	var
		dwords, i, ubn, ubw:Int32;
		j, last, laststep, q, t, stp, carry:Int32;
		xd, xn, rund:double;
		resul: array of double;
		n: array of double;
		d: array of double;
		w: array of double;
		fac1, fac2:decfloat;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;
		ubn := 2 * dwords + 3+14;
		ubw := ubn + 4+1;
		fac1 := x;
		fac2 := y;

		if (fac2.exponent = 0) then // if fac2 = 0, return
		begin
			// a divide-by-zero error and
			// bail out.
			for i := 0 to dwords do
				fac1.mantissa[i] := 99999999;

			fac1.exponent := 99999 + bias + 1;
			//er=DIVZ_ERR;
			exit(fac1);
		end;
		if (fac1.exponent = 0) then //fact1=0, just return
			exit(fac1);
		setLength(w, ubw + 1);
		setLength(resul, ubn + 5);
		setLength(n, ubn + 5);
		setLength(d, ubn + 5);
		
		//double w[ubw];

		for j := 0 to dwords do
		begin
			n[2 * j + 2] := fac1.mantissa[j] div 10000;
			n[2 * j + 3] := fac1.mantissa[j] mod 10000;
			d[2 * j + 2] := fac2.mantissa[j] div 10000;
			d[2 * j + 3] := fac2.mantissa[j] mod 10000;
		end;
		n[1] := (fac1.exponent and $7fffffff) - bias - 1;
		d[1] := (fac2.exponent and $7fffffff) - bias - 1;
		for j := ubn to ubw do
			w[j] := 0;
		t := ubn - 1;
		w[1] := n[1] - d[1] + 1;
		w[2] := 0;
		for j := 2 to ubn do
			w[j + 1] := n[j];
		xd := (d[2] * b + d[3]) * b + d[4] + d[5] / b;
		laststep := t + 2;
		for stp := 1 to laststep do
		begin
			xn := ((w[(stp + 2) - 1] * 10000 + w[(stp + 2)]) * 10000 + w[(stp + 2) + 1]) * 10000;
			if (ubw >= ((stp + 2) + 2)) then
				xn := xn + w[(stp + 2) + 2];
			q := round(xn / xd);
			last := min(stp + t + 1, ubw);
			//subtract(w, q, d, (stp + 2), last);
			for j := (stp + 2) to last do
				w[j] := w[j] - q * d[j - (stp + 2) + 2];
			//normalize(w, (stp + 2), q);
			w[(stp + 2)] := w[(stp + 2)] + w[(stp + 2) - 1] * 10000;
			w[(stp + 2) - 1] := q;
		
		end;
		//finalnorm(w, (laststep + 1));
		for j := (laststep + 1) downto 3 do
		begin
			if w[j] < 0 then
				carry := round((-w[j] - 1) / b) + 1
			else
			begin
				if w[j] >= b then
					carry := -round(w[j] / b)
				else
					carry := 0;
			end;
			w[j] := w[j] + carry * b;
			w[j - 1] := w[j - 1] - carry;
		end;
		
		
		
		
		If (w[2] <> 0) Then
			dec(laststep);
		rund := w[laststep + 1] / b;
		w[laststep] := w[laststep];
		If (rund >= 0.5) Then
			w[laststep]+=1;
		if (w[2] = 0) then
		begin
			for j := 1 to (t + 1) do
				resul[j] := w[j + 1];
		end
		else
		begin
			for j := 1 to (t + 1) do
				resul[j] := w[j];
		end;
// left-shift out leading 0's if any to preserve accuracy
{if 1}

		j := 0;
		if (resul[2] < 10) then
			j := 3
		else if (resul[2] < 100) then
			j := 2
		else if (resul[2] < 1000) then
			j := 1;

		if (j > 0) then
			LSHIFT_da(resul, 2 * dwords + 4, j);
{endif}
		resul[1] := w[1];
		If w[2] = 0 Then
			resul[1]-=1;

		for j := 0 to dwords do
			fac1.mantissa[j] := round(resul[2 * j + 2] * 10000 + resul[2 * j + 3]);

		j := norm_fac1(fac1, dwords);
		fac1.exponent := round(resul[1] + bias);

		fac1.sign := fac1.sign xor fac2.sign;
		result := fac1;
	end;

	function fpdiv_si(const num:decfloat; const den:Int32; const dwords_in:Int32 = num_dwords) : decfloat;
	var
		dwords, i:Int32;
		fac1, fac2:decfloat;
		carry, remder:uint64;
		divisor, quotient:int64;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;
		remder:=0;
		if (den < 0) then
			divisor := -den
		else
			divisor := den;

		fac1 := num;
		if (divisor = 0) then
		begin
			//printf('%s\n', 'error: divisor = 0');
			for i := 0 to dwords do result.mantissa[i] := 99999999;
			result.exponent := 99999999;
			result.sign := num.sign;
			exit(result);
		end;
		if (divisor > 99999999) then
		begin
			fac2 := si2fp(den, dwords);
			result := fpdiv(fac1, fac2, dwords);
			exit(result);
		end;

		for i := 0 to (dwords) do
		begin
			quotient := fac1.mantissa[i] + remder * 100000000;
			remder := quotient mod divisor;
			fac1.mantissa[i] := quotient div divisor;
		end;
		quotient := remder * 100000000;
		quotient := quotient div divisor;
		carry := fac1.mantissa[0];

		if (carry = 0) then
		begin
			lshift_8(fac1, dwords);
			fac1.exponent -= 8;
			fac1.mantissa[dwords] += Int32(quotient)
		end
		else if (carry < 10) then
		begin
			lshift_7(fac1, dwords);
			fac1.exponent -= 7;
			fac1.mantissa[dwords] += Int32(quotient div 10);
		end
		else if (carry < 100) then
		begin
			lshift_6(fac1, dwords);
			fac1.exponent -= 6;
			fac1.mantissa[dwords] += Int32(quotient div 100);
		end
		else if (carry < 1000) then
		begin
			lshift_5(fac1, dwords);
			fac1.exponent -= 5;
			fac1.mantissa[dwords] += Int32(quotient div 1000);
		end
		else if (carry < 10000) then
		begin
			lshift_4(fac1, dwords);
			fac1.exponent -= 4;
			fac1.mantissa[dwords] += Int32(quotient div 10000);
		end
		else if (carry < 100000) then
		begin
			lshift_3(fac1, dwords);
			fac1.exponent -= 3;
			fac1.mantissa[dwords] += Int32(quotient div 100000);
		end
		else if (carry < 1000000) then
		begin
			lshift_2(fac1, dwords);
			fac1.exponent -= 2;
			fac1.mantissa[dwords] += Int32(quotient div 1000000);
		end
		else if (carry < 10000000) then
		begin
			lshift_1(fac1, dwords);
			fac1.exponent -= 1;
			fac1.mantissa[dwords] += Int32(quotient div 10000000);
		end;

		i := norm_fac1(fac1, dwords);
		if (den < 0) then fac1.sign := fac1.sign xor (-1);

		//i=norm_fac1(fac1, dwords);
		result := fac1;
	end;

	function fpfrac(const num:decfloat; const dwords_in:Int32 = num_dwords) : decfloat;
	var
		dwords : Int32;
		n : decfloat;

	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;
		n:=fpfix( num, dwords);
		result := fpsub(num, n, dwords);
	end;

	//returns the positive of n
	Function fpabs (const n:decfloat): decfloat;
	var
		x: decfloat;
	begin
		x:=n;
		x.sign:=0;
		Result := x;
	End;

	//changes the sign of n, if n is positive then n will be negative & vice versa
	Function fpneg (const n: decfloat): decfloat;
	var
		x: decfloat;
	begin
		x:=n;
		x.sign:=x.sign Xor -1;
		Result := x;
	End;

	function fp2str_fix (const n:decfloat; const places_in : Int32=num_digits) : ansiString;
	var
		places, i, ex : Int32;
		v, ts, s :ansistring;
	begin
		places:=places_in;
		If (places>(num_digits)) Then places:=num_digits;

		if n.exponent > 0 then
			ex := (n.exponent and $7fffffff) - bias - 1
		else
			ex := 0;

		if n.sign<>0 then
			s := '-'
		else
			s := ' ';

		ts := trim(IntToStr(n.mantissa[0]));
		if length(ts) < 8 then ts := ts + stringn(8 - length(ts), '0');
		v := ts;
		for i := 1 to num_dwords do
		begin
			ts := trim(IntToStr(n.mantissa[i]));
			if length(ts) < 8 then ts := stringn(8 - length(ts), '0') + ts;
			v := v + ts;
		end;
		if places < num_digits then v := leftstr(v, places);
		if ex = 0 then
			v := leftstr(v, 1) + '.' + rightstr(v, length(v)-1)
		else if ex < 0 then
			v := '.' + stringn(abs(ex) - 1, '0') + v
		else if ex > 0 then
			v := leftstr(v, ex + 1) + '.' + rightstr(v, length(v)-1-ex);

		result := (s + v).trimright(['0']);
		result := result.trimright(['.']);
	end;

	function fp2str(const n:decfloat; const places_in : Int32=90) : ansiString;
	var
		dwords, ex, e, ee, sign, places : Int32;
		z:decfloat;
	begin
		places:=places_in;
		if places>num_digits-1 then
			places:=num_digits-2
		else
			places+=1;
		dwords:=num_dwords;
		if n.exponent <> 0 then
			ex := (n.exponent and $7fffffff) - bias - 1
		else
			ex := 0;

		if n.exponent=0 then exit(' 0');
		e:=ex;
		ee:=0;
		z:=si2fp (5, dwords);
		sign:=n.sign;
		if (ex > (-5)) and (e < (places-2)) then
		begin
			if e<0 then ee:=1;
			z.exponent := (-(places-e-ee) + bias + 1);
			if sign=0 then
				z:=fpadd(z, n, dwords)
			else
				z:=fpsub(n, z, dwords);

			exit(fp2str_fix(z, places-ee));
		end
		else
		begin
			z.exponent := (-(places-e-1) + bias + 1);
			if sign=0 then
				z:=fpadd(z, n, dwords)
			else
				z:=fpsub(n, z, dwords);

			exit(fp2str_exp(z, places-1));
		end;
	end;

	function fp2str2(const n:decfloat; const places_in : Int32=num_digits) : ansiString;
	var
		s : ansiString;
	begin
		s:=fp2str(n, places_in);
		if length(s)<places_in then
			s:=s+stringn(places_in-length(s), ' ');
		Result := s;
	end;

	function ipower(const x:double; const e:Int32) : double;
	var
		z, y : double;
		n : Int32;
	begin
		//take x to an integer power
		y := x;
		n := abs(e);
		z := 1.0;
		while (n > 0) do
		begin
			while ((n and 1) = 0) do
			begin
				n := n div 2;
				y := y * y;
			end;
			n := n - 1;
			z := y * z;
		end;
		if (e < 0) then z := 1.0 / z;
		result := z;
	end;

	Function log2_64(value : Uint64) : int32;
	//https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
	const tab64 : array [0..63] of int32 = (63, 0, 58, 1, 59, 47, 53, 2, 60, 39, 48, 27, 54, 33, 42,
		3, 61, 51, 37, 40, 49, 18, 28, 20, 55, 30, 34, 11, 43, 14, 22, 4, 62, 57, 46, 52, 38, 26, 32,
		41, 50, 36, 17, 19, 29, 10, 13, 21, 56, 45, 25, 31, 35, 16, 9, 12, 44, 24, 15, 8, 23, 7, 6, 5);
	begin
		value := value Or Uint64(value Shr 1);
		value := value Or Uint64(value Shr 2);
		value := value Or Uint64(value Shr 4);
		value := value Or Uint64(value Shr 8);
		value := value Or Uint64(value Shr 16);
		value := value Or Uint64(value Shr 32);
		log2_64 := tab64[(Uint64((value - (value Shr 1)) * $07EDD5E59A4E28C2) Shr 58)];
	End;

	function fpsqr(const num:decfloat; const dwords_in:Int32 = num_dwords) : decfloat;
	var
		dwords, ex, k, l, prec:Int32;
		r, r2, tmp, n:decfloat;
		x, y:double;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		//l := round(log2((dwords + 1) * 8 * double(0.0625))) + 1;
		l := Int32(round(Ln(double(((num_digits+8)*0.0625)*1.5))));
		
		//l=estimated number of iterations needed
		//first estimate is accurate to about 16 digits
		//l is approximatly = to log2((NUM_DIGITS+9)/16)
		//NUM_DIGITS+9 because decfloat has an extra 9 guard digits
		n := num;
		tmp:=si2fp(0, dwords);

		if ((fpcmp(n, tmp, dwords) = 0) or (fpcmp(n, tmp, dwords) < 0)) then
		begin
			result:=si2fp(0, dwords);
			exit(result)
		end;
		result:=si2fp(1, dwords);
		if (fpcmp(n, result, dwords) = 0) then exit(result);

		//=====================================================================
		//hack to bypass the limitation of double exponent range
		//in case the number is larger than what a double can handle
		//for example, if the number is 2e500
		//we separate the exponent and mantissa in this case 2
		//if the exponent is odd then multiply the mantissa by 10
		//take the square root and assign it to decfloat
		//divide the exponent in half for square root
		//in this case 1.414213562373095e250
		
		if (n.exponent > 0) then
			ex := (n.exponent and $7fffffff) - bias - 1
		else
			ex := 0;

		x := double(n.mantissa[0]) + (double(n.mantissa[1]) / double(100000000));
		x := x / double(10000000);

		if (x = 0) then
		begin
			result:=si2fp(0, dwords);
			exit(result);
		end;
		if ((abs(ex) mod 2) <> 0) then
		begin
			x := x * 10;
			ex -= 1;
		end;
		x := sqrt(x); //approximation
		y := x * 10000000;
		r:=si2fp(0, dwords);
		r.mantissa[0] := round(y);
		r.mantissa[1] := round((y - r.mantissa[0]) * 100000000);
		r.exponent := ex div 2 + bias + 1;
		r.sign:=0;
		
		//if len(v)>1 and k=0 then r.exponent+=1

		//=====================================================================
		//newton-raphson method
		prec := 3;
		for k := 1 to (l + 2) do
		begin
			prec := 2 * prec - 1;
			tmp:=fpdiv(n, r, prec);
			r2:=fpadd(r, tmp, prec);
			r:=fpdiv_si(r2, 2, prec);
		end;
		
		result := r;
	end;

	procedure pi_chudnovsky_binary_split(const a:uint64; const b:uint64; var pab:decfloat; var qab:decfloat; var tb:decfloat; const dwords_in:Int32 = num_dwords);
	var
		dwords : Int32;
		c545140134, c13591409, temp, c3_over_24 : decfloat;
		pam, qam, tam, pmb, qmb, tmb:decfloat;
		m, t:uint64;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		c545140134:=si2fp( 545140134, dwords);
		c13591409:=si2fp(13591409, dwords);
		c3_over_24:=si2fp(0, dwords);
		c3_over_24.exponent := 1073741841;
		c3_over_24.mantissa[0] := 10939058;
		c3_over_24.mantissa[1] := 86003200;

		if ((b - a) = 1) then
		begin
			if (a = 0) then
			begin
				pab:=si2fp(1, dwords);
				qab:=si2fp(1, dwords);
			end
			else
			begin
				t := 6 * a;
				pab:=si2fp(t - 1, dwords);
				pab:=fpmul_si( pab, a + a - 1, dwords);
				pab:=fpmul_si(pab, t - 5, dwords);
				temp:=si2fp(a, dwords);
				qab:=fpmul(temp, temp, dwords);
				qab:=fpmul(qab, temp, dwords);
				qab:=fpmul(qab, c3_over_24, dwords);
			end;
			copydf(tb, c545140134, dwords);
			tb:=fpmul_si(tb, a, dwords);
			tb:=fpadd(tb, c13591409, dwords);
			tb:=fpmul(tb, pab, dwords);
			if ((a mod 2) <> 0) then
				tb.sign := (tb.sign xor (-1));
		end
		else
		begin
			m := (a + b) div 2;
			pi_chudnovsky_binary_split(a, m, pam, qam, tam, dwords);
			pi_chudnovsky_binary_split(m, b, pmb, qmb, tmb, dwords);
			pab:=fpmul(pam, pmb, dwords);
			qab:=fpmul(qam, qmb, dwords);
			temp:=fpmul(qmb, tam, dwords);
			tb:=fpmul(pam, tmb, dwords);
			tb:=fpadd(tb, temp, dwords);
		end;
	end;

	function pi_chudnovsky_bs(const digits_in:Int32 = num_digits) : decfloat;
	var
		dwords, digits : Int32;
		p, q, sqrtc, t : decfloat;
		n:uint64;
		digits_per_term:double;
	begin
		dwords := digits_in div 8;
		if (dwords > num_dwords) then dwords := num_dwords;
		digits:=8*dwords;
		digits_per_term := 14.18164746272548;
		n := round(digits / digits_per_term + 1);
		pi_chudnovsky_binary_split(0, n, p, q, t, dwords);
		sqrtc:=si2fp(10005, dwords);
		sqrtc:=fpsqr( sqrtc, dwords);
		result:=fpmul_si(sqrtc, 426880, dwords);
		result:=fpmul(result, q, dwords);
		result:=fpdiv(result, t, dwords);
	end;

	function fpipow(const x:decfloat; const e:int64; const dwords_in:Int32 = num_dwords) : decfloat;
	var
		dwords : Int32;
		one, y, r:decfloat;
		n:int64;

	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		//take x to an long power

		copydf(y, x, dwords);
		if (e < 0) then
			n := -e
		else
			n := e;

		r:=si2fp(1, dwords);
		copydf(one, r, dwords);
		while (n > 0) do
		begin
			while ((n and 1) = 0) do
				begin
				n:=n shr 1;
				y:=fpmul(y, y, dwords);
				//c += 1;
			end;
			n -= 1;
			r:=fpmul(y, r, dwords);
			//c += 1;
		end;
		if (e < 0) then	r:=fpdiv(one, r, dwords);
		result:=r;
	end;

    procedure MyFrExp(
      x: Double;
      var y: Double;
      var e: Int32
    );

    var
      xUInt: UInt64;
      ExpBin: Int32;

    begin
      xUInt := pUInt64( @x )^;             // transfer Double to UInt64;
      ExpBin := ( xUint shr 52 ) and $7ff; // get biased exponent

      if( ( ExpBin<>0 ) and ( ExpBin<>$7ff ) ) then begin // sort out zero, de-normal, NaN and Inf cases
        xUInt := ( xUInt and UInt64( $800fffffffffffff ) ) or $3fe0000000000000;
        y := pDouble( @xUInt )^;
        e := ExpBin - $3fe;
      end else begin
        if( ( ExpBin=0 ) and ( ( xUInt and $fffffffffffff )=0 ) ) then begin
          e := 0;
          y := pDouble( @xUInt )^;
        end else begin
          // none of the remaining values - de-normal, NaN or Inf - can begin
          // 'normalised' in a meaningfull way. Add whatever you like in error
          // handling here
        end;
      end;
    end;

	function dbl2fp(const x:double; const dwords_in:Int32 = num_dwords) : decfloat;
	var
		dwords, n : Int32;
		z:double;
		y:extended;
		fp, w:decfloat;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;
		result:=si2fp(0, dwords);
		Frexp(abs(x), y, n);
		z := y * 100000000;
		result.mantissa[0] := trunc(z);
		z := z - result.mantissa[0];
		z := z * 100000000;
		result.mantissa[1] := trunc(z);
		result.exponent := bias;
		result.sign := 0;
		norm_fac1(result, 2);
		if (x < 0) then result.sign := -1;
		fp:=si2fp(2, dwords);
		w:=fpipow(fp, n, dwords);
		result:=fpmul(result, w, dwords);
	end;

	function fp2i(const x:decfloat) : Int32;
	var
		e, ex, n:Int32;
		tmp:decfloat;
	begin
		tmp:=si2fp(2147483647, 2);
		if (fpcmp_abs(x, tmp, 3)>0) then exit(2147483647);
		if (x.exponent <> 0 ) then
			ex:=(x.exponent and $7fffffff)-bias-1
		else
			ex:=0;

		e:=abs(ex);
		tmp:=si2fp(5, 3);
		tmp.exponent-=1;
		tmp:=fpadd(tmp, x, 3);
		tmp:=fpfix(tmp, 2);
		if (e>7) then
			n:=tmp.mantissa[0]*Int32(p10_constants[e-7])+(tmp.mantissa[1] div Int32(p10_constants[15-e]))
		else
			n:=tmp.mantissa[0] div Int32(p10_constants[7-e]);

		if (x.sign<>0) then n:=-n;
		result := n;
	end;

	function fp2ui(const x:decfloat) : longword;
	var
		e, ex:Int32;
		m:longword;
		tmp:decfloat;
	begin
		tmp:=si2fp(0);
		if (fpcmp(x, tmp, 3)<0) then exit(0);
		tmp:=si2fp(4294967295, 2);
		if (fpcmp(x, tmp, 3)>0) then exit(4294967295);
		if (x.exponent <> 0 ) then
			ex:=(x.exponent and $7fffffff)-bias-1
		else
			ex:=0;

		e:=abs(ex);
		tmp:=si2fp(5, 3);
		tmp.exponent-=1;
		tmp:=fpadd(tmp, x, 3);
		tmp:=fpfix(tmp, 2);
		if (e>7) then
			m:=tmp.mantissa[0]*longword(p10_constants[e-7])+(tmp.mantissa[1] div longword(p10_constants[15-e]))
		else
			m:=tmp.mantissa[0] div longword(p10_constants[7-e]);

		result := m;
	end;

	function fp2ui64(const x:decfloat) : uint64;
	var
		e, ex:Int32;
		m:uint64;
		tmp:decfloat;
	begin
		tmp:=si2fp(0);
		if (fpcmp(x, tmp, 3)<0) then exit(0);

		tmp.mantissa[0]:=18446744;
		tmp.mantissa[1]:=7370955;
		{$if num_dwords>1} tmp.mantissa[2]:=16150000; {$endif}
		tmp.exponent:=19+bias+1;

		if (fpcmp(x, tmp, 3)>0) then exit(18446744073709551615);
		if (x.exponent <> 0 ) then
			ex:=(x.exponent and $7fffffff)-bias-1
		else
			ex:=0;

		e:=abs(ex);
		tmp:=si2fp(5, 3);
		tmp.exponent-=1;
		tmp:=fpadd(tmp, x, 3);
		tmp:=fpfix(tmp, 3);

		if (e>15) then
		begin
			m:=(uint64(tmp.mantissa[0])*100000000+uint64(tmp.mantissa[1]))*uint64(p10_constants[e-15]);
			{$if num_dwords>1} m:=m+(uint64(tmp.mantissa[2]) div uint64(p10_constants[23-e])); {$endif}
		end
		else if (e>7) then
			m:=uint64(tmp.mantissa[0])*uint64(p10_constants[e-7])+uint64(tmp.mantissa[1]) div uint64(p10_constants[15-e])
		else
			m:=uint64(tmp.mantissa[0]) div uint64(p10_constants[7-e]);
		result := m;
	end;

	function ui2fp(const m : uint64; const dwords_in : Int32=num_dwords) : decfloat;
	var
		dwords, ex, i : Int32;
		mm : uint64;
	begin
		dwords:=dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;
		mm:=m;
		for i := 0 to dwords do result.mantissa[i] := 0;

		if (mm = 0) then
		begin
			result.exponent := 0;
			result.sign := 0;
			exit(result);
		end;
		ex := Int32(log10_64(mm));
		result.exponent := bias + ex + 1;

		if (ex > 15) then
		begin
			ex -= 7;
			result.mantissa[0] := (mm div uint64(p10_constants[ex]));
			mm := mm mod uint64(p10_constants[ex]);
			ex -= 8;
			result.mantissa[1] := (mm div uint64(p10_constants[ex]));
			mm := mm mod uint64(p10_constants[ex]);
			{$if num_dwords>1} result.mantissa[2] := (mm * uint64(p10_constants[8 - ex]));{$endif}
		end
		else if (ex >= 8) then
		begin
			ex -= 7;
			result.mantissa[0] := (mm div uint64(p10_constants[ex]));
			mm := mm mod uint64(p10_constants[ex]);
			result.mantissa[1] := (mm * uint64(p10_constants[8 - ex]));
		end
		else if (ex < 8) then
			result.mantissa[0] := (mm * uint64(p10_constants[7 - ex]));
	end;


//=====================================================================\\
	function GetTickCount : DWORD;
	{On Windows, this is number of milliseconds since Windows was
	started. On non-Windows platforms, LCL returns number of
	milliseconds since Dec. 30, 1899, wrapped by size of DWORD.
	This value can overflow Int32 variable when checks turned on,
	so "wrap" value here so it fits within Int32.
	Also, since same thing could happen with Windows that has been
	running for at least approx. 25 days, override it too.}
	begin
		{$IFDEF MSWINDOWS}
			Result := Windows.GetTickCount mod High(Int32);
		{$ELSE}
			Result := LclIntf.GetTickCount mod High(Int32);
		{$ENDIF}
	end;


//OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\\
	procedure writelnq(const n:decfloat; const digits : Int32=56);
	begin
		writeln(fp2str_exp(n, digits));
	end;

	procedure writeq(const n:decfloat; const digits : Int32=56);
	begin
		write(fp2str_exp(n, digits));
	end;

	function decfloat.ToString(places:Int32 = num_digits): ansistring;
	begin
		Result:=fp2str(Self, places);
	end;

	function decfloat.ToStringExp(places:Int32 = num_digits): ansistring;
	begin
		Result:=fp2str_exp(Self, places);
	end;

	function decfloatc.ToString(places:Int32 = num_digits): ansistring;
	begin
		Result:=fp2str(Self.re, places);
		if Self.im < 0 then
			Result:=Result+' -'+fp2str(abs(Self.im), places)
		else
			Result:=Result+' +'+fp2str(Self.im, places);
	end;

	function decfloatc.ToStringExp(places:Int32 = num_digits): ansistring;
	begin
		Result:=fp2str_exp(Self.re, places);
		if Self.im < 0 then
			Result:=Result+' -'+fp2str_exp(abs(Self.im), places)
		else
			Result:=Result+' +'+fp2str_exp(Self.im, places);
	end;

	operator := (r : ansistring) z : decfloat;

	begin
	  z:=str2fp(r);
	end;

	operator := (r : decfloat) z : ansistring;

	begin
	  z:=fp2str(r);
	end;

	operator := (r : Int16) z : decfloat;

	begin
	  z:=si2fp(Int32(r));
	end;

	operator := (r : Int32) z : decfloat;

	begin
	  z:=si2fp(r);
	end;

	operator := (r : decfloat) z : Int32;

	begin
	  z:=fp2i(r);
	end;

	operator := (r : int64) z : decfloat;

	begin
	  z:=si2fp(r);
	end;

	operator := (r : double) z : decfloat;

	begin
	  z:=dbl2fp(r);
	end;

	operator + (x : decfloat; y : decfloat) z : decfloat;

	begin
	  z:=fpadd(x, y);
	end;

	operator - (x : decfloat) z : decfloat;

	begin
	  z:=x;
	  z.sign := (z.sign xor (-1));
	end;

	operator - (x : decfloat; y : decfloat) z : decfloat;

	begin
	  z:=fpsub(x, y);
	end;

	operator * (x : decfloat; y : decfloat) z : decfloat;

	begin
	  z:=fpmul(x, y);
	end;

	operator * (x : decfloat; y : Int32) z : decfloat;

	begin
	  z:=fpmul_si(x, y);
	end;

	operator / (x : decfloat; y : decfloat) z : decfloat;

	begin
	  z:=fpdiv(x, y);
	end;

	operator / (x : decfloat; y : Int32) z : decfloat;

	begin
	  z:=fpdiv_si(x, y);
	end;

	Operator ** (const x : decfloat; const y : decfloat) : decfloat;
	begin
		result := exp(log(x)*y);
	End;

	Operator ** (const x : Double; const y : decfloat) : decfloat;
	begin
		result := exp(Log(dbl2fp(x))*y);
	End;

	Operator ** (const x : Integer; const y : decfloat) : decfloat;
	begin
		result := exp(Log(si2fp(x))*y);
	End;
	
	operator < (x : decfloat; y : decfloat) z : boolean;
	var c:Int32;
	begin
	  c:=fpcmp(x, y);
	  z:=c<0;
	end;

	operator <= (x : decfloat; y : decfloat) z : boolean;
	var c:Int32;
	begin
	  c:=fpcmp(x, y);
	  z:=c<=0;
	end;

	operator > (x : decfloat; y : decfloat) z : boolean;
	var c:Int32;
	begin
	  c:=fpcmp(x, y);
	  z:=c>0;
	end;

	operator >= (x : decfloat; y : decfloat) z : boolean;
	var c:Int32;
	begin
	  c:=fpcmp(x, y);
	  z:=c>=0;
	end;

	operator <> (x : decfloat; y : decfloat) z : boolean;
	var c:Int32;
	begin
	  c:=fpcmp(x, y);
	  z:=c<>0;
	end;

	operator = (x : decfloat; y : decfloat) z : boolean;
	var c:Int32;
	begin
	  c:=fpcmp(x, y);
	  z:=c=0;
	end;

	function fpsign( const x : decfloat): Int32;
	begin
		fpsign:=0;
		if x<0 then fpsign:=-1;
		if x=0 then fpsign:=0;
		if x>0 then fpsign:=1;
	end;

	function FloatToStr(x : decfloat):string; overload;
	begin
		FloatToStr:=fp2str(x);
	end;

	procedure Str(x : decfloat; var s:string); overload;
	begin
		s:=fp2str(x);
	end;

	function sqrt(x : decfloat): decfloat;overload;
	begin
		result:=fpsqr(x);
	end;

	function abs(x : decfloat): decfloat;overload;
	begin
		result:=x;
		result.sign:=0;
	end;

	function sin(x : decfloat): decfloat;overload;
	begin
		result:=fpSin(x);
	end;

	function cos(x : decfloat): decfloat;overload;
	begin
		result:=fpCos(x);
	end;

	function tan(x : decfloat): decfloat;overload;
	begin
		result:=fpTan(x);
	end;

	function arcsin(x : decfloat): decfloat;overload;
	begin
		result:=fpArcSin(x);
	end;

	function arccos(x : decfloat): decfloat;overload;
	begin
		result:=fpArcCos(x);
	end;

	function arctan(x : decfloat): decfloat;overload;
	begin
		result:=fpArcTan(x);
	end;

	function log(x : decfloat): decfloat;overload;
	begin
		result:=log_agm(x);
	end;

	function exp(x : decfloat): decfloat;overload;
	begin
		result:=fpexp(x);
	end;

	function frac(x : decfloat): decfloat;overload;
	begin
		frac:=fpfrac(x);
	end;

	function trunc(x : decfloat): decfloat;overload;
	begin
		trunc:=fpfix(x);
	end;
//=======================================================================================

	function fpSin(x: decfloat): decfloat;

	const	A0  = '1';
			A1  = '-.166666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666656';
			A2  = '0.833333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333312837003e-2';
			A3  = '-0.198412698412698412698412698412698412698412698412698412698412698412698412698412698412698412698412697737444966048e-3';
			A4  = '0.275573192239858906525573192239858906525573192239858906525573192239858906525573192239858906525484671365838005109e-5';
			A5  = '-0.250521083854417187750521083854417187750521083854417187750521083854417187750521083854417181588192427454858764465e-7';
			A6  = '0.160590438368216145993923771701549479327257105034882812660590438368216145993923771701285849508268402610988318809e-9';
			A7  = '-0.764716373181981647590113198578807044415510023975632441240906849372457838066303599107532637130446982142423935374e-12';
			A8  = '0.281145725434552076319894558301032001623349273520453103397392224033991852228715216952979475805635178847750668098e-14';
			A9  = '-0.822063524662432971695598123687228074922073899182611413442667321736818048347828018574106418544315177713042338812e-17';
			A10 = '0.195729410633912612308475743735054303552874737900621765105396981363228345228295916277591726315408380188118367037e-19';
			A11 = '-0.386817017063068403771691193152281232317934264625734713647029367346232525194225170575352200281115325527816415253e-22';
			A12 = '0.644695028438447339619485321920468720529890441042891189394707950182291117405447388211089657694464302812356450942e-25';
			A13 = '-0.918368986379554614842571683647391339786168719434316969735598055943905677408894019425312585423397367812757768506e-28';
			A14 = '0.113099628864477169315587645769383169924405014704255241809975058402524330880458392137874957162344724446324394808e-30';
			A15 = '-0.121612504155351794962997468569229214972478347403208386145390651489209475069212788822882821279979140347180989975e-33';
			A16 = '0.115163356207719502805868814932982211143284272893708437522400970268434523750106762605826106348737771802376046945e-36';
			A17 = '-0.967759295863189099208981638092286297524381942605107581375324863884120569249596678608433115725091393475246935864e-40';
			A18 = '0.726546017915307131538274503048930928819093512253493206563834798135763419653380326778941958110535786470998875789e-43';
			A19 = '-0.490246975651354339769415628184136528258723207005956239555439070006973503563879662185349639566666696134572393516e-46';
			A20 = '0.298931082714240451074404483096711329413689860034259328570612315176298156298993813519951441431266463635873707713e-49';
			A21 = '-0.165521086774219476155120850023448611714109095583805285809579158799143086443557827034675155030521545645379394110e-52';
			A22 = '0.835965084715278702393674259072136165997911382223366608479978580515691863814526765694836292192450284180870442060e-56';
			A23 = '-0.386662836669782516462665275221862042794934626023349938316999004328812238872208204848547336560505964670320988852e-59';
			A24 = '0.164352468623969149943788481677979108436615041101655682622684231056139699136404729500407839377392414007591137747e-62';

	var sn,x2: decfloat;
		sign:Int32;
	begin
		sign:=0;
		if x<0 then
		begin
			sign:=-1;
			x:=-x;
		end;
		if x>fpPi2 then
		begin
			x:=fpfrac(x/fpPi2);
			x:=x*fpPi2;
		end;
		if x>fpPihalf
			then
				if x<Pi3o2
					then
						x:=fpPi-x
					else
						x:=x-fpPi2;
		x:=x/9;
		x2:=x*x;
		sn:=(A0+(A1+(A2+(A3+(A4+(A5+(A6+(A7+(A8+(A9+(A10+(A11+(A12+(A13+(A14+(A15+
			(A16+(A17+(A18+(A19+(A20+(A21+(A22+(x2*A24+A23)*x2)*x2)*x2)*x2)*x2)*x2)
			*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x;

	  if sign<0 then sn:=-sn;
	  sn:=sn*(3-4*sn*sn);

	  fpSin:=sn*(3-4*sn*sn)
	end;

	function fpCos(x: decfloat): decfloat;

	const	A0  = '1';
			A1  = '-.4999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999994721';
			A2  = '0.4166666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666665621360871e-1';
			A3  = '-0.1388888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888854451270574738e-2';
			A4  = '0.2480158730158730158730158730158730158730158730158730158730158730158730158730158730158730158725644213389781219540e-4';
			A5  = '-2.755731922398589065255731922398589065255731922398589065255731922398589065255731922398588750980614595596512318923e-7';
			A6  = '2.087675698786809897921009032120143231254342365453476564587675698786809897921009032106698274232523330624155210153e-9';
			A7  = '-1.147074559772972471385169797868210566623265035963448661861360274058686757099455126283818763190613529714553377891e-11';
			A8  = '4.779477332387385297438207491117544027596937649847702757755667808577861487835681124988711847586663858811113753943e-14';
			A9  = '-1.561920696858622646221636435005733342351940408446961685541067911299953546188895284977828153146115130065583351444e-16';
			A10 = '4.110317623312164858477990618436140374610369495913057067213336608547373892481324660795350325932672520623662928130e-19';
			A11 = '-8.896791392450573286748897442502468343312488086391898413881668726420645861793068520597680120778186378293173759925e-22';
			A12 = '1.611737571096118349048713304801171801324726102607227973442549983675870905814278441599245567654898745663177847333e-24';
			A13 = '-2.479596263224797460074943545847956617422655542472653504987407428935336567238877943440604966558279635032505664489e-27';
			A14 = '3.279889237069837910152041727312111927807745426326503439096280743521705083405323141198625400715500889291751485579e-30';
			A15 = '-3.769987628815905643852921525646105664146825508872167896410951976022887555989319206462275922154423847904269702672e-33';
			A16 = '3.800390754854743592593670892788412967640253420811360978105269101130113351836803478383682018285354576347489926649e-36';
			A17 = '-3.387157535521161847231435733323000135617276359939640722898889509770724867519203645280731260776540092274152525090e-39';
			A18 = '2.688220266286636386691615661248345947271627617428900493013208786573066818885796824441820180098748129679804128498e-42';
			A19 = '-1.911963205040281925100720510987268731741680241998275528414009592120512524012702432579358988121591591326830083814e-45';
			A20 = '1.225617439128385849400550899069151514579508607057473027440126990112219354212293319416478760765428830433984001016e-48';
			A21 = '-7.117406731291437132972785958274464124745589390878139936192587433963538109320745624028616476864082913696927070763e-52';
			A22 = '3.761842881216953274249787675345772704624526037429344479318689496767530861512138234511295249094362011904925941139e-55';
			A23 = '-1.817315326458027900185306127402970583114969176062051137517511979800387650439115823948958513219179624857155393161e-58';
			A24 = '8.053180970768418501403321884065945547834282740190199439315571732533516060616256383923068265144713840919403240939e-62';

	var sn,x2: decfloat;
		sign:Int32;
	begin
		sign:=1;
		if x<0 then
		begin
			x:=-x;
		end;
		if x>fpPi2 then
		begin
			x:=fpfrac(x/fpPi2);
			x:=x*fpPi2;
		end;
		if x>fpPihalf
			then
				if x<Pi3o2
					then
					begin
						x:=fpPi-x;
						sign:=-sign;
					end
					else
						x:=x-fpPi2;

		x:=x/9;
		x2:=x*x;
		sn:=A0+(A1+(A2+(A3+(A4+(A5+(A6+(A7+(A8+(A9+(A10+(A11+(A12+(A13+(A14+(A15+
			(A16+(A17+(A18+(A19+(A20+(A21+(A22+(x2*A24+A23)*x2)*x2)*x2)*x2)*x2)*x2)
			*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2;

	  if sign<0 then sn:=-sn;
	  sn:=sn*(4*sn*sn-3);
	  sn:=sn*(4*sn*sn-3);
	  fpCos:=sn;
	end;

	function fpTan(x: decfloat): decfloat;

	const	A0 =' 1';
			A1='.3381385686116249656720414612193072552971221238659188023840633055582528630195389012407063167533882134699597287926106479';
			A2='-0.3047840831546141695108148800286182936222558831019111115275865277367853433854498337060905525973361168111386564380498399e-1';
			A3='0.8912143674122393190443198653257534845486948027622749934280165068266358307592609669379625164314141385255120006985075964e-3';
			A4='-0.1210241064430942027435700958617037890048906053736323329070677276326957628880530460187540982611238638995983442435319156e-4';
			A5='8.874462637694657633897988585961093227825235906348754864224172682040758465728601535879443567358601898913435194514293259e-8';
			A6='-3.789554665160215407670561703948534243666431556104788762769288702493249641531182352063993383208477581346800389612056156e-10';
			A7='9.758052094374283728664660196424204830820381310518517364607177916759732410655218783605454239561211542692708585893711999e-13';
			A8='-1.523031736882553015227321360324491385056548499874296923304544577506761082925793825617263359159972945143935402758203148e-15';
			A9='1.406025150227028131991569012973774250653659056316990654868116504244059294581111028544718876852895169047553715763583716e-18';
			A10='-7.195368052693418126232336227900639214993773760639840332158698857192405158415061476783993100587828347680069791341622604e-22';
			A11='1.775196338612049895421224428292158941865542822591701483498991849866864483457750499506820085377545047955426806615823160e-25';
			A12='-1.492797878342368213063833730815962654744748080789666364681698643306003635735762985340399731094878158557636350899725174e-29';
			A13='1.047656960440587288899711204892311824544206008757357333496066353530501199811093243004418846066487111616356461148741256e-36';

			B1='1.014415705834874897016124383657921765891366371597756407152189916674758589058616703722118950260164640409879186377830135';
			B2='-.4972015072803342096596942174717541944432233135696758963191519249909390386390816316006747458832666912072932714437766914';
			B3='0.3731551268872398992316193684958684280090466395842519377854264253141841805270582636702643205095521521308469514253286869e-1';
			B4='-0.1017478601466026385713768005220505119975023836834532403464341680297768129689218413915816823395242239653511376576719367e-2';
			B5='0.1331048998319742210236800464449450117954129064404712037693059254150717144044422285829187637081732197452515013053642962e-4';
			B6='-9.533150933975377628157227276535992804054345360928033244720319207497643291334208107228286391827264442227852744345548253e-8';
			B7='4.004388942947812626075143801207961189956156207849522998398555945631870533159506408235101037719029444068891028000856936e-10';
			B8='-1.018451721752856124433949066919449718990812399072155073031609379408521594394859755518896250474474171971312405316404202e-12';
			B9='1.574030468133099942010478564846302168780742667781844056959794954997650782373689287027812086770267730267167680877972800e-15';
			B10='-1.441214300339117113226436108185397989879439645494110323387795564326960702693412679972728925513123618658075614684594501e-18';
			B11='7.322683225767239269064513888188599547196653123038598273263877801423691825514541452550034128321487855766497060401049634e-22';
			B12='-1.794732773821548650740289175272993412159426422485019922713166384787682484963298596588267224885478276801865322279779300e-25';
			B13='1.499232567538698718668462157852134215984693421744738381345540432220872216574724975454687329166622654604277992650738348e-29';

	var tn,x2: decfloat;
		sign :Int32;
	begin
		sign:=1;
		if x<0 then
		begin
			sign:=-1;
			x:=-x;
		end;

		if x>fpPi2 then
		begin
			x:=fpfrac(x/fpPi2);
			x:=x*fpPi2;
		end;

		while x >  fpPihalf do x := x - fpPi;
		while x <= -fpPihalf do x := x + fpPi;

		x:=x/9;
		x2:=x*x;
		tn:=x+x2*x*(A1+(A2+(A3+(A4+(A5+(A6+(A7+(A8+(A9+(A10+(A11+(A12+A13*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)/
		(B1+(B2+(B3+(B4+(B5+(B6+(B7+(B8+(B9+(B10+(B11+(B12+B13*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2)*x2);

		if sign<0 then tn:=-tn;

		tn:=(tn*(3-tn*tn))/(A0-3*tn*tn);
		tn:=(tn*(3-tn*tn))/(A0-3*tn*tn);
		fpTan:=tn;
	end;

	function fpArcTan(x_in: decfloat): decfloat;

	var at, k1, k2, x2, c, d, h, m, s, x: decfloat;
			e, f, g, i, n, s1, s2, sg: integer;

	begin
	  x:=x_in;
	  s1:=fpSign(x);
	  x:=Abs(x);
	  if x<1 then
		  begin
			s2:=1;
			k1:=0
		  end
		else
		  begin
			x:=1/x;
			s2:=-1;
			k1:=fpPihalf
		  end;
	  if x>fpTMS3
		then
		  begin
			x:=(x*fpSQR3-1)/(x+fpSQR3);
			k2:=fpPIo6
		  end
		else
		  k2:=0;
	  if x = 0 then
		   fpArcTan:=0
		 else
		   begin
			 if Abs(x) < 4e-5 then
			   n:=1
			 else
			   n:=fp2i(fpfix((480/7)*Abs(x)+150/7));
			 x:=1/x;
			 x2:=x*x;
			 c:=0;
			 s:=0;
			 d:=(4*n+1)*x-1/x;
			 m:=2*(x-1/x);
			 f:=4*n*n;
			 g:=(2*n-1);
			 e:=4*g;
			 h:=1;
			 for i:=1 to 2*n-1 do
			   h:=h*x;
			 for i:=1 to n do
			   begin
				 s:=1/(g*h)-s;
				 c:=f/(c+d);
				 f:=f-e;
				 g:=g-2;
				 e:=e-8;
				 h:=h/x2;
				 d:=d-m
			   end;
			 c:=c+d;
			 for i:=1 to n do
				 c:=c*x2;
			 sg:=Trunc(1-2* (n Mod 2));
			 at:=(s+sg/c);
			 fpArcTan:=(k1+s2*(k2+at))*s1
		   end;
	end;

	function fpArcSin(x: decfloat): decfloat;

	begin
	  fpArcSin:=0;
	  if Abs(x)<=1
		then
		  if Abs(x)<1
			then
			  fpArcSin:=fpArcTan(x/(Sqrt(1-x*x)))
			else
			  fpArcSin:=fpSign(x)*fpPihalf
		else
		  begin
			Write(^G);
			WriteLn('fpArcSin Error')
		  end;
	end;

	function fpArcCos(x: decfloat): decfloat;

	begin
	  fpArcCos:=0;
	  if Abs(x)<=1
		then
		  if Abs(x)<1
			then
			  fpArcCos:=2*fpArcTan(Sqrt((1-x)/(1+x)))
			else
			  fpArcCos:=fpPihalf*(1-x)
		 else
		   begin
			 Write(^G);
			 WriteLn('fpArcCos Error')
		   end;
	end;

//===================================== complex =========================================
// the code is a translation of the C++ code published on the web with some tweaks to make the results agree with Mathematica
// the file has since vanished from the web but his revised work is available at https://github.com/RobTillaart/Complex/tree/master
// 
//     FILE: Complex.cpp
//   AUTHOR: Rob Tillaart
//  VERSION: 0.1.09
//  PURPOSE: library for Complex math for Arduino
//      URL: http://arduino.cc/playground/Main/ComplexMath
// 
//  Released to the public domain

(*
	operator := (r : ansistring) z : decfloatc;

	begin
	  z:=str2fp(r);
	end;
*)
	function Str2c(const s: ansistring): decfloatc;
	var
		a, b: string;
		c, i, k, l: LongInt;
		Sgn1 (*, Sgn2 *): LongInt;
		resultC: decfloatc;
	begin
		a := UpperCase(Trim(s));
		Sgn1 := 1;
		//Sgn2 := 1;
		resultC.re := 0;
		resultC.im := 0;

		if (Length(a) = 0) or (a[Length(a)] <> 'I') then
		begin
			WriteLn('not a valid complex number');
			Exit(resultC);
		end;

		if (a = 'I') or (a = '+I') then
		begin
			resultC.im := 1;
			Exit(resultC);
		end
		else if (a = '-I') then
		begin
			resultC.im := -1;
			Exit(resultC);
		end;

		l := Length(a);
		i := l - 1;              { Pascal strings are 1-based }

		while (i >= 1) and (a[i] = ' ') do
			Dec(i);

		b := a[i];

		if b <> '*' then
		begin
			if (b = '+') or (b = '-') then
			begin
				if b = '-' then
					resultC.im := -1
				else
					resultC.im := 1;
				Inc(i);
			end
			else
			begin
				WriteLn('not a valid complex number');
				Exit(resultC);
			end;
		end;

		a := Trim(Copy(a, 1, i - 1));

		if (Length(a) > 0) and (a[1] = '-') then
		begin
			Sgn1 := -1;
			Delete(a, 1, 1);
		end;

		l := Length(a);
		c := Pos('E', a);

		if c > 0 then
		begin
			k := c + 1;
			while (k <= l) and not (a[k] in ['+', '-']) do
				Inc(k);
		end
		else
		begin
			k := 1;
			while (k <= l) and not (a[k] in ['+', '-']) do
				Inc(k);
		end;

		//if (k > 1) and (a[k - 1] = '-') then
		//	Sgn2 := -1;

		resultC.re := Sgn1 * Str2fp(Copy(a, 1, k - 1));

		if resultC.im = 0 then
		resultC.im := Str2fp(Copy(a, k, MaxInt));

		Result := resultC;
	end;

	operator := (r : decfloatc) z : ansistring;

	begin
		z:=fp2str(r.re, 105);
		if r.im <0 then
			z:=z+' -'+fp2str(abs(r.im), 105)+'*I'
		else
			z:=z+' +'+fp2str(r.im, 105)+'*I';
	end;
(*
	operator := (r : Int16) z : decfloatc;

	begin
		z.re:=si2fp(Int32(r));
	end;
*)
(*
	operator := (r : Int32) z : decfloatc;

	begin
		z.re:=si2fp(r);
	end;
*)
(*
	operator := (r : decfloatc) z : Int32;

	begin
		z:=fp2i(r.re);
	end;
*)
(*
	operator := (r : int64) z : decfloatc;

	begin
		z.re:=si2fp(r);
	end;
*)
(*
	operator := (r : double) z : decfloatc;

	begin
		z.re:=dbl2fp(r);
	end;
*)
//***********************************************
	Operator + (x, y : decfloatc) z : decfloatc;
	begin
		z.re := x.re + y.re;
		z.im := x.im + y.im;
	End;

	Operator + (x : decfloatc; y : decfloat) z : decfloatc;
	begin
		z.re := x.re + y;
		z.im := x.im;
	End;

	Operator + (x : decfloat; y : decfloatc) z : decfloatc;
	begin
		z.re := x + y.re;
		z.im := y.im;
	End;

	Operator - (y : decfloatc) z : decfloatc; //negate
	begin
		z.re := -y.re;
		z.im := -y.im;
	End;

	Operator - (x, y : decfloatc) z : decfloatc;
	begin
		z.re := x.re - y.re;
		z.im := x.im - y.im;
	End;

	Operator - (x : decfloatc; y : Double) z : decfloatc;
	begin
		z.re := x.re - y;
		z.im := x.im;
	End;

	Operator - (x : Double; y : decfloatc) z : decfloatc;
	begin
		z.re := x - y.re;
		z.im := -y.im;
	End;

	Operator * (x, y : decfloatc) z : decfloatc;
	begin
		z.re := x.re * y.re - x.im * y.im;
		z.im := x.re * y.im + x.im * y.re;
	End;

	Operator * (x : decfloatc; y : decfloat) z : decfloatc;
	begin
		z.re := x.re * y;
		z.im := x.im * y;
	End;

	Operator * (x : decfloat; y : decfloatc) z : decfloatc;
	begin
		z.re := x * y.re;
		z.im := x * y.im;
	End;


	Operator / (x, y : decfloatc) z : decfloatc;
	var
		f : decfloat;
	begin
		f := 1/(y.re*y.re + y.im*y.im);
		z.re := (x.re * y.re + x.im * y.im) * f;
		z.im := (x.im * y.re - x.re * y.im) * f;
	End;

	Operator / (x : decfloatc; y : Double) z : decfloatc;
	var f : decfloat;
	begin
		f := 1/(y*y);
		z.re := (x.re * y) * f;
		z.im := (x.im * y) * f;
	End;

	Operator / (x : Double; y : decfloatc) z : decfloatc;
	var f : decfloat;
	begin
		f := 1/(y.re*y.re + y.im*y.im);
		z.re := (x * y.re) * f;
		z.im := (-x * y.im) * f;
	End;

	Function csquare(x : decfloatc) : decfloatc;
	begin
		Result.re := x.re * x.re - x.im * x.im;
		Result.im := 2 * x.re * x.im;
	End;

	Function reciprocal(x :decfloatc) : decfloatc;
	var
		f, r, i :decfloat;
	begin
		f := 1/(x.re*x.re + x.im*x.im);
		r := x.re*f;
		i := -x.im*f;
		Result.re := r;
		Result.im := i;
	End;

	Function cabs(x : decfloatc) : decfloatc;
	begin
		Result.re := Sqrt(x.re*x.re+x.im*x.im);
		Result.im :=0;
	End;

	function ArcTan2(const Y, X: decfloat): decfloat; overload;
	begin
		if (X = 0) then
		begin
			if (Y = 0) then Result := 0
			else if (Y > 0) then Result := fpPihalf
			else Result := -fpPihalf;
		end
		else
		begin
			Result := ArcTan(Y / X);
			if (X < 0) then
			begin
				if (Y < 0) then Result := Result - fppi
				else Result := Result + fppi;
			end;
		end;
	end;

	function cArctan2(y, x: decfloatc): decfloatc; overload;
	var
		num1, den, ratio: decfloatc;
		log_result: decfloatc;
	begin
		// arctan2(y, x) = -i * log((x + i*y) / sqrt(x^2 + y^2))

		// Calculate numerator: x + i*y
		num1.re := x.re - y.im;	// x.re + i*y.im => real: x.re, imag: y.im
		num1.im := x.im + y.re;	// x.im + i*y.re => real: y.re, imag: x.im

		// Calculate denominator: sqrt(x^2 + y^2)
		den.re := x.re * x.re - x.im * x.im + y.re * y.re - y.im * y.im;	// real part of x^2 + y^2
		den.im := 2 * x.re * x.im + 2 * y.re * y.im;	// imag part of x^2 + y^2
		den := cSqrt(den);

		// Avoid division by zero
		if (abs(den.re) < '1e-105') and (abs(den.im) < '1e-105') then
		begin
			// Special case: denominator is zero
			if (abs(num1.re) < '1e-105') and (abs(num1.im) < '1e-105') then
			begin
				// Both numerator and denominator are zero - undefined
				Result.re := 0; //NaN;
				Result.im := 0; //NaN;
				Exit;
			end
			else
			begin
				// Only denominator is zero
				Result.re := 0; // NaN;
				Result.im := 0; //NaN;
				Exit;
			end;
		end;

		// Calculate ratio: (x + i*y) / sqrt(x^2 + y^2)
		// Complex division: (a+bi)/(c+di) = ((ac+bd)/(c^2+d^2)) + ((bc-ad)/(c^2+d^2))i
		ratio.re := (num1.re * den.re + num1.im * den.im) / (den.re * den.re + den.im * den.im);
		ratio.im := (num1.im * den.re - num1.re * den.im) / (den.re * den.re + den.im * den.im);

		// Take complex logarithm
		log_result := cLog(ratio);

		// (-i) * log(ratio) = -i*(a+bi) = -i*a - i^2*b = b - i*a
		Result.re := log_result.im;	// b
		Result.im := -log_result.re; // -a
	end;

	function hypot(const Y, X: decfloat): decfloat; overload;
	begin
		Result := (Sqrt(((x)*(x))+((y)*(y))));
	end;

	function modulus(const Y, X: decfloat): decfloat; overload; // same as hypot
	begin
		Result := (Sqrt(((x)*(x))+((y)*(y))));
	end;

	function phase(const Y, X: decfloat): decfloat; overload;
	begin
		Result := ArcTan2(x,y);
	end;

	Function csqrt(x : decfloatc) : decfloatc;
	var m : decfloat;
	begin
		m := hypot(x.re,x.im);
		Result.re := Sqrt(0.5 * (m+x.re));
		Result.im := Sqrt(0.5 * (m-x.re));
		If x.im < 0 Then Result.im := -Result.im;
	End;

	Function cexp(x : decfloatc) : decfloatc;
	var e : decfloat;
	begin
		e := Exp(x.re);
		Result.re := e * Cos(x.im);
		Result.im := e * Sin(x.im);
	End;

	Function clog(x : decfloatc) : decfloatc;
	var m, p : decfloat;
	begin
		m := modulus(x.re, x.im);
		p := phase(x.re, x.im);
		If p > fppi Then p := p-fppi2;
		Result.re := Log(m);
		result.im := p;
	End;

	Function cpow(Const x, c : decfloatc) : decfloatc;
	begin
		Result := cexp(clog(x)*c);
	End;

	Function clogn(Const x, c : decfloatc) : decfloatc;
	begin
		Result := clog(x)/clog(c);
	End;

	Function clog10(Const x : decfloatc) : decfloatc;
	begin
		Result := clog(x)*fpconst1olog10;
	End;

	Function sinh(const x : decfloat) : decfloat; overload;
	var y : decfloat;
	begin
		y := Exp(x);
		Result := (y-1/y)*0.5;
	End;

	Function cosh(const x : decfloat) : decfloat; overload;
	var y : decfloat;
	begin
		y := Exp(x);
		Result := (y+1/y)*0.5;
	End;

	Function tanh(const x : decfloat) : decfloat; overload;
	var y : decfloat;
	begin
		y := Exp(x);
		Result := (y-1/y)/(y+1/y);
	End;

	Function csin(const z : decfloatc) : decfloatc;
	var i : decfloatc;
	begin
		i.re := 0; i.im := 1;
		Result := -i/2*(cexp(i*z)-cexp(-i*z));
	End;

	Function ccos(const z : decfloatc) : decfloatc;
	var i : decfloatc;
	begin
		i.re := 0; i.im := 1;
		Result := (cexp(i*z)+cexp(-i*z))/2;
	End;

	Function ctan(const z : decfloatc) : decfloatc;
	var i : decfloatc;
	begin
		i.re := 0; i.im := 1;
		Result := (cexp(i*z) - cexp(-i*z)) / (i*(cexp(i*z) + cexp(-i*z)));
	End;

	Function gonioHelper1(const x : decfloatc; const mode : Byte) : decfloatc;
	var c, i : decfloatc;
	begin
		i.re := 0; i.im := -1;
		c := csqrt(fponec - csquare(x));
		If mode = 0 Then
			c := c + x * i
		Else
			c := x + c * i;
		i.im := 1;
		Result := clog(c) * i;
	End;

	Function cArcSin(const x : decfloatc) : decfloatc;
	begin
		Result := gonioHelper1(x, 0);
	End;

	Function cArcCos(const x : decfloatc) : decfloatc;
	begin
		Result := gonioHelper1(x, 1);
	End;

	Function cArcTan(const x : decfloatc) : decfloatc;
	var mi, a, b : decfloatc;
	begin
		mi.re := 0; mi.im := -1;
		a.re := x.re; a.im := x.im - 1;
		b.re := -x.re; b.im := -x.im - 1;
		Result := (mi * clog(a/b)) * 0.5;
	End;

	Function ccsc(const x : decfloatc) : decfloatc;
	begin
		Result := fponec / csin(x);
	End;

	Function csec(const x : decfloatc) : decfloatc;
	begin
		Result := fponec / ccos(x);
	End;

	Function ccot(const x : decfloatc) : decfloatc;
	begin
		Result := fponec / ctan(x);
	End;

	Function cArcCsc(const x : decfloatc) : decfloatc;
	begin
		Result := fponec / cArcSin(x);
	End;

	Function cArcSec(const x : decfloatc) : decfloatc;
	begin
		Result := fponec / cArcCos(x);
	End;

	Function cArcCot(const x : decfloatc) : decfloatc;
	begin
		Result := fponec / cArcTan(x)
	End;

	//HYPERBOLICUS I
	Function csinh(const x : decfloatc) : decfloatc;
	var y : decfloatc;
	begin
		y := cexp(x);
		Result := (y-reciprocal(y))*0.5;
	End;

	Function ccosh(const x : decfloatc) : decfloatc;
	var y : decfloatc;
	begin
		y := cexp(x);
		Result := (y+reciprocal(y))*0.5;
	End;

	Function ctanh(const x : decfloatc) : decfloatc;
	begin
		Result := csinh(x) / ccosh(x);
	End;

	Function gonioHelper2(const x : decfloatc; const mode : Byte) : decfloatc;
	var c : decfloatc;
	begin
		c := csquare(x);
		If mode = 0 Then
			c := c+1
		Else
			c := c-1;

		Result := clog(x + csqrt(c));
	End;

	Function cArcSinh(const x : decfloatc) : decfloatc;
	begin
		Result := gonioHelper2(x,0);
	End;

	Function cArcCosh(const x : decfloatc) : decfloatc;
	begin
		Result := clog(x+csqrt(x-fponec)*csqrt(x+fponec));
		//return gonioHelper2(x,1)
	End;

	Function cArcTanh(const x : decfloatc) : decfloatc;
	var c : decfloatc;
	begin
		c := clog(x + fponec);
		c := c - clog(-(x - fponec));
		Result := c * 0.5;
	End;

	//HYPERBOLICUS II
	Function ccsch(const x : decfloatc) : decfloatc;
	begin
		Result := fponec / csinh(x);
	End;

	Function csech(const x : decfloatc) : decfloatc;
	begin
		Result := fponec / ccosh(x);
	End;

	Function ccoth(const x : decfloatc) : decfloatc;
	begin
		Result := fponec / ctanh(x);
	End;

	Function cArcCsch(const x : decfloatc) : decfloatc;
	begin
		Result := cArcSinh(fponec/x);
	End;

	Function cArcSech(const x : decfloatc) : decfloatc;
	begin
		Result := cArcCosh(fponec/x);
	End;

	Function cArcCoth(const x : decfloatc) : decfloatc;
	begin
		Result := cArcTanh(fponec/x);
	End;

	Function fpsignc(const x : decfloatc) : Int32;
	var
		sign : Int32;
	begin
		sign := fpSign(x.re);
		If sign=0 Then
			sign := fpSign(x.im);

		Result :=sign;
	End;

	Operator ** (const x : decfloatc; const y : decfloatc) : decfloatc;
	begin
		result := cexp(clog(x)*y);
	End;

	Operator ** (const x : Double; const y : decfloatc) : decfloatc;
	begin
		result := cexp(Log(dbl2fp(x))*y);
	End;

	Operator ** (const x : Integer; const y : decfloatc) : decfloatc;
	begin
		result := cexp(Log(si2fp(x))*y);
	End;

	function toDecfloatc(const x : decfloat; const y : decfloat) : decfloatc; overload;
	begin
		Result.re :=x;
		Result.im :=y;
	end;

	function toDecfloatc(const x : double; const y : decfloat) : decfloatc; overload;
	begin
		Result.re :=x;
		Result.im :=y;
	end;

	function toDecfloatc(const x : decfloat; const y : double) : decfloatc; overload;
	begin
		Result.re :=x;
		Result.im :=y;
	end;

	function toDecfloatc(const x : double; const y : double) : decfloatc; overload;
	begin
		Result.re :=x;
		Result.im :=y;
	end;

//
	function toDecfloatc(const x : Int32; const y : decfloat) : decfloatc; overload;
	begin
		Result.re :=x;
		Result.im :=y;
	end;

	function toDecfloatc(const x : decfloat; const y : Int32) : decfloatc; overload;
	begin
		Result.re :=x;
		Result.im :=y;
	end;

	function toDecfloatc(const x : Int32; const y : Int32) : decfloatc; overload;
	begin
		Result.re :=x;
		Result.im :=y;
	end;


function fpgamma(x : decfloat) : decfloat;
var z : decfloat;
	k, n, m : int32;
begin
	result:=0;
	if (x=1) then
	begin
		result:=1;
		exit;
	end;

	if (x<0) then
	begin
		n:=fp2i(x);
		if (x<>n) then
		begin
			result:=fppi/(fpgamma(1-x)*sin(fppi*x));
		end
		else
		begin
			result:='1e99999999';
		end;
	end
	else
	begin
		n:=fp2i(abs(x));
		m:=0;
		if n<169 then
			m:=169-n;
		x:=x+m;
		z:=1/x;
		result:=x*(gamma_coeff[0]+(gamma_coeff[1]+(gamma_coeff[2]+(gamma_coeff[3]+(gamma_coeff[4]+
		(gamma_coeff[5]+(gamma_coeff[6]+(gamma_coeff[7]+(gamma_coeff[8]+(gamma_coeff[9]+
		(gamma_coeff[10]+(gamma_coeff[11]+(gamma_coeff[12]+(gamma_coeff[13]+(gamma_coeff[14]+
		(gamma_coeff[15]+(gamma_coeff[16]+(gamma_coeff[17]+(gamma_coeff[18]+(gamma_coeff[19]+
		(gamma_coeff[20]+(gamma_coeff[21]+(gamma_coeff[22]+(gamma_coeff[23]+(gamma_coeff[24]+
		(gamma_coeff[25]+(gamma_coeff[26]+(gamma_coeff[27]+(gamma_coeff[28]+(gamma_coeff[29]+
		(gamma_coeff[30]+(gamma_coeff[31]+(gamma_coeff[32]+(gamma_coeff[33]+(gamma_coeff[34]+
		(gamma_coeff[35]+(gamma_coeff[36]+(gamma_coeff[37]+(gamma_coeff[38]+(gamma_coeff[39]+
		(gamma_coeff[40]+(gamma_coeff[41]+(gamma_coeff[42]+(gamma_coeff[43]+(gamma_coeff[44]+
		(gamma_coeff[45]+(gamma_coeff[46]+(gamma_coeff[47]+(gamma_coeff[48]+(gamma_coeff[49]+
		(gamma_coeff[50]+(gamma_coeff[51]+(gamma_coeff[52]+(gamma_coeff[53]+(gamma_coeff[54]+
		(gamma_coeff[55]+(gamma_coeff[56]+(gamma_coeff[57]+(gamma_coeff[58]+(gamma_coeff[59]+
		(gamma_coeff[60]+(gamma_coeff[61]+(gamma_coeff[62]+(gamma_coeff[63]+(gamma_coeff[64]+
		(gamma_coeff[65]+(gamma_coeff[66]+(gamma_coeff[67]+(gamma_coeff[68]+(z*gamma_coeff[70]+
		gamma_coeff[69])*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*
		z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*
		z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z);
		
		result:=sqrt2pi*z*sqrt(z)/(z**(x)*exp(x))*result;
		
		if m>0 then
		begin
			for k:=1 to m do
			begin
				z:=x-k;
				result:=result/z;
			end;
		end;
	end;
end;

function fpgammac(x : decfloatc) : decfloatc;
var z : decfloatc;
	y : decfloat;
	k, n, m : int32;
begin
	result.re:=0; result.im:=0;

	if (x.re=0) and (x.im=0) then
	begin
		result.re:=1;
		exit;
	end;

	if (x.re<0) then
	begin
		n:=fp2i(x.re);
		if (x.re<>n) or (x.im<>0) then
		begin	
			z.re:=fppi; z.im:=0;
			result:=z/(fpgammac(1-x)*csin(fppi*x));
		end
		else
		begin
			result.re:='1e99999999';
			result.im:='1e99999999';
		end;
	end
	else
	begin

		y:=sqrt(x.re * x.re + x.im * x.im);
		n:=fp2i(y);
		m:=0;
		if n<169 then
			m:=169-n;
		x:=x+m;
		z:=1/x;
		result:=x*(gamma_coeff[0]+(gamma_coeff[1]+(gamma_coeff[2]+(gamma_coeff[3]+(gamma_coeff[4]+
		(gamma_coeff[5]+(gamma_coeff[6]+(gamma_coeff[7]+(gamma_coeff[8]+(gamma_coeff[9]+
		(gamma_coeff[10]+(gamma_coeff[11]+(gamma_coeff[12]+(gamma_coeff[13]+(gamma_coeff[14]+
		(gamma_coeff[15]+(gamma_coeff[16]+(gamma_coeff[17]+(gamma_coeff[18]+(gamma_coeff[19]+
		(gamma_coeff[20]+(gamma_coeff[21]+(gamma_coeff[22]+(gamma_coeff[23]+(gamma_coeff[24]+
		(gamma_coeff[25]+(gamma_coeff[26]+(gamma_coeff[27]+(gamma_coeff[28]+(gamma_coeff[29]+
		(gamma_coeff[30]+(gamma_coeff[31]+(gamma_coeff[32]+(gamma_coeff[33]+(gamma_coeff[34]+
		(gamma_coeff[35]+(gamma_coeff[36]+(gamma_coeff[37]+(gamma_coeff[38]+(gamma_coeff[39]+
		(gamma_coeff[40]+(gamma_coeff[41]+(gamma_coeff[42]+(gamma_coeff[43]+(gamma_coeff[44]+
		(gamma_coeff[45]+(gamma_coeff[46]+(gamma_coeff[47]+(gamma_coeff[48]+(gamma_coeff[49]+
		(gamma_coeff[50]+(gamma_coeff[51]+(gamma_coeff[52]+(gamma_coeff[53]+(gamma_coeff[54]+
		(gamma_coeff[55]+(gamma_coeff[56]+(gamma_coeff[57]+(gamma_coeff[58]+(gamma_coeff[59]+
		(gamma_coeff[60]+(gamma_coeff[61]+(gamma_coeff[62]+(gamma_coeff[63]+(gamma_coeff[64]+
		(gamma_coeff[65]+(gamma_coeff[66]+(gamma_coeff[67]+(gamma_coeff[68]+(z*gamma_coeff[70]+
		gamma_coeff[69])*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*
		z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*
		z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z)*z);
		
		result:=sqrt2pi*z*csqrt(z)/(z**x*cexp(x))*result;
	
		if m>0 then
		begin
			for k:=1 to m do
			begin
				z:=x-k;
				result:=result/z;
			end;
		end;
	end;

end;

Function fpzetac(n : decfloatc) : decfloatc;
	//Riemann zeta function evaluated using the Euler-Maclaurin summation formula.

	//first we sum the first 70 terms of the Zeta series,
	//then use the Euler-Maclaurin summation formula
var	
	nc, m, i, k, k1, sign : Int32;
	s, dx, n1, one, two, half, ic : decfloatc;
begin
	sign:=0;
	nc:=70;
	result.re:=Int32(0);
	result.im:=Int32(0);
	one.re:=Int32(1);
	one.im:=Int32(0);
	two.re:=Int32(2);
	two.im:=Int32(0);
	half.re:=dbl2fp(0.5);
	half.im:=Int32(0);
	If n.re<0 Then
	begin
		n:=1-n;
		sign:=-1;
		if (fpfrac(n.re)=0) and (n.im=0) then
		begin
			m:=fp2i(abs(n.re));
			if odd(m) then
				exit(result); // complex(0,0)
		end;
	End;
	
	//  f(k)=(k+nc)^-n

	//  inf           nc              inf                        inf
	//  ====          ====           /                          ==== B
	//  \       1     \       1      [                          \    (2*k)   (2*k-1)
	//   >    ----- =  >    ----- +  I    f(k) dk + B * f(0) -   >   ------ f  (0)
	//  /      n      /      n       ]               1          /    (2*k)!
	//  ====  k       ====  k       /                           ====
	//  k = 1         k = 1          0                          k = 1

	s.re:=Int32(0);
	s.im:=Int32(0);
	
	For i:=2 To nc do
	begin
		ic.re:=i;
		ic.im:=0;
		s:=s+ic**(-n);
	end;
	
	ic.re:=nc;
	ic.im:=0;
	//                            inf
	//                          /
	//                         [
	s:=s+ic**(one-n)/(n-one); //      I    f(k) dk  = nc^(1-n)/(n-1)
	//                         ]
	//                         /
	//                          0

	dx:=-n*ic**(-n-one);    //first derivative of f(0)
	s:=s-half*ic**(-n);     //f(0)*B(1), first term of the Euler-Maclaurin summation formula.
	n1:=one/(ic*ic);
	k:=1;
	k1:=2;
	For i:=1 To 75 do
	begin
		s:=s-dx*Bernoulli_fac[i];
	  
			 //next derivative
							// (2*i-1)
		dx:=dx*n1*(n+k1)*(n+k); //f  (0)
		k:=k+2;
		k1:=k1+2;
	end;
	s:=one+s;
	ic.re:=fpPi2;
	ic.im:=Int32(0);
	If sign<0 Then
		s:=two/ic**n*ccos(n*fpPihalf)*fpgammac(n)*s;

	Result := s;
End;

Function fpzeta(n : decfloat) : decfloat;
	//Riemann zeta function evaluated using the Euler-Maclaurin summation formula.

	//first we sum the first 70 terms of the Zeta series,
	//then use the Euler-Maclaurin summation formula
var	
	nc, m, i, k, k1, sign : Int32;
	s, dx, n1, one, two, ic : decfloat;
begin
	sign:=0;
	nc:=70;
	result:=Int32(0);
	one:=Int32(1);
	two:=Int32(2);

	If n<0 Then
	begin
		n:=1-n;
		sign:=-1;
		if (fpfrac(n)=0) and (n=0) then
		begin
			m:=fp2i(abs(n));
			if odd(m) then
				exit(result);
		end;
	End;
	
	//  f(k)=(k+nc)^-n

	//  inf           nc              inf                        inf
	//  ====          ====           /                          ==== B
	//  \       1     \       1      [                          \    (2*k)   (2*k-1)
	//   >    ----- =  >    ----- +  I    f(k) dk + B * f(0) -   >   ------ f  (0)
	//  /      n      /      n       ]               1          /    (2*k)!
	//  ====  k       ====  k       /                           ====
	//  k = 1         k = 1          0                          k = 1

	s:=Int32(0);
	
	For i:=2 To nc do
	begin
		ic:=i;
		s:=s+ic**(-n);
	end;
	
	ic:=nc;
	
	//                            inf
	//                          /
	//                         [
	s:=s+ic**(one-n)/(n-one); //      I    f(k) dk  = nc^(1-n)/(n-1)
	//                         ]
	//                         /
	//                          0

	dx:=-n*ic**(-n-one);    //first derivative of f(0)
	s:=s-0.5*ic**(-n);     //f(0)*B(1), first term of the Euler-Maclaurin summation formula.
	n1:=one/(ic*ic);
	k:=1;
	k1:=2;
	For i:=1 To 75 do
	begin
		s:=s-dx*Bernoulli_fac[i];
	  
			 //next derivative
							// (2*i-1)
		dx:=dx*n1*(n+k1)*(n+k); //f  (0)
		k:=k+2;
		k1:=k1+2;
	end;
	
	s:=one+s;

	If sign<0 Then
		s:=two/fpPi2**n*fpcos(n*fpPihalf)*fpgamma(n)*s;

	Result := s;
End;

//===================================== end complex =====================================

	function agm(const a_in:decfloat; const b_in:decfloat; const digits:Int32=num_digits) : decfloat;
	var
		a, b, c, d, eps:decfloat;
		s:string;
	begin
		a:=a_in;
		b:=b_in;
		str(digits-5, s);
		s:=trim(s);
		eps:='1e-'+s;
		repeat
			c := (a + b)/2;
			d := sqrt(a*b);
			a := c;
			b := d;
		until abs(a-b)<eps;
		result:=a;
	end;

	function calcln2(const digits:Int32=num_digits) : decfloat;
	begin
		result:=decfloat(2);
		result:=fpipow(result, -2*digits+1);
		result:=(fppi/(2*agm(1, result, digits)))/(2 * digits + 1);
	end;

	function log_agm(const x:decfloat; const digits:Int32=num_digits) : decfloat;
	var
		y:decfloat;
	begin

		if x <= 0 then
		begin
			writeln('can nott take fplog_agm of zero or a negative number');
			exit(0);
		end;

		if x=1 then exit(0);

		if x < 1 then
		begin
			y:=2;
			y:=fpipow(y, 2 - 2 * digits);
			y*=x;
			y:=fppi/(2*agm(1, y, digits))-fpln2 * 2 * digits;
			y.sign:=-1;
		end
		else
		begin
			y:=2;
			y:=fpipow(y, 2 - 2 * digits);
			y/=x;
			y:=fppi/(2*agm(1, y, digits))-fpln2 * 2 * digits;
		end;
		result:=y
	end;

	function fpexp(const x2:decfloat; const dwords_in:Int32 = num_dwords):decfloat;
	const	A0 ='1';
			A1 ='1';
			A2 ='.5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000658819';
			A3 ='.1666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666524824002543';
			A4 ='0.4166666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666717502466219295548e-1';
			A5 ='0.8333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333326071623994270402980903e-2';
			A6 ='0.1388888888888888888888888888888888888888888888888888888888888888888888888888888888888888888944127087045005952179329605e-2';
			A7 ='0.1984126984126984126984126984126984126984126984126984126984126984126984126984126984126981534497883198069903565331823523e-3';
			A8 ='0.2480158730158730158730158730158730158730158730158730158730158730158730158730158730240783874862051041215243115881770734e-4';
			A9 ='0.2755731922398589065255731922398589065255731922398589065255731922398589065255730064216498435903309349274761934797635460e-5';
			A10='2.755731922398589065255731922398589065255731922398589065255731922398589065287127437691754620532423297643861581554400507e-7';
			A11='2.505210838544171877505210838544171877505210838544171877505210838544171469343342041026059445393714486361837569624727690e-8';
			A12='2.087675698786809897921009032120143231254342365453476564587675698790989427448401798711874753882229604378151083784860402e-9';
			A13='1.605904383682161459939237717015494793272571050348828126605904349360501600396942306488084607334943877364547242551497514e-10';
			A14='1.147074559772972471385169797868210566623265035963448661861589480044898141895699826060866995077269484032965287156898938e-11';
			A15='7.647163731819816475901131985788070444155100239756324399825345751274554759233560323654386311742879321099647185407412010e-13';
			A16='4.779477332387385297438207491117544027596937649847760026779373598006483516127715098440006711848035345240759140895842950e-14';
			A17='2.811457254345520763198945583010320016233492734987149493427524600124974167468815526772681746615806969279779090586917115e-15';
			A18='1.561920696858622646221636435005733342351941099566558357437710055927460580547362009137909023332293200586700492981466003e-16';
			A19='8.220635246624329716955981236872280749202289832546264473848011704065438699178721726401133664247102105613303398929641775e-18';
			A20='4.110317623312164858477990618436140415991660445867089735532136485552727343983082918905034346332847547651043702628442535e-19';
			A21='1.957294106339126123084757437350465130349568713230362628098826339876429930203677445710161578664447307741019733573689543e-20';
			A22='8.896791392450573286748897443729731301067008869437849768321366311016482355949846725861892370253131047024434755143106281e-22';
			A23='3.868170170630684037716910322380374135100948034490266570735919485786248145040164377043511889088531707500761915322081371e-23';
			A24='1.611737571096118349050455032684633927859537257199385729326861492510874762285658367187972751493057541445554834997361396e-24';
			A25='6.446950284384473380816184523176777688136691115175486294514818141695854291969051735719770439225898685711320583675775938e-26';
			A26='2.479596263224808348503459452771326010737844016105437398546695413856302900407937506262708691304643261209264762990232689e-27';
			A27='9.183689863735254024933144217132060143387792480293931734798962193831600268278075142181163556662643757563028833596425494e-29';
			A28='3.279889262207765240241404732834255550489877284696521209106124544712712920926169874727225671702016117308250957706385696e-30';
			A29='1.130988870786209636795984578198190902818273592522597443112638497729739938583224740666276313994130608289748440058308742e-31';
			A30='3.783789750286940938141798620984947851907247605475528683305672908635712394061165443956859922287453079262950148925726142e-33';


	var
		dwords, i, n:Int32;
		fac, x, temp, accum:decfloat;
	begin
		dwords := dwords_in;
		if (dwords > num_dwords) then dwords := num_dwords;

		temp:=si2fp(0, dwords);
		fac:=si2fp(1, dwords);

		if (fpcmp(x2, temp, dwords) = 0) then
		begin
			copydf(result, fac, dwords);
			exit(result);
		end;
		fac:=fpfix(abs(x2), dwords);
		x:=abs(x2)-fac;
		if x>0 then
		begin
			x:=fpdiv_si(x, 128, dwords);
			accum:=A0+(A1+(A2+(A3+(A4+(A5+(A6+(A7+(A8+(A9+(A10+(A11+(A12+(A13+(A14+(A15+(A16+
				(A17+(A18+(A19+(A20+(A21+(A22+(A23+(A24+(A25+(A26+(A27+(A28+(A29+A30*x)*x)*x)*x)
				*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x;

			for i:=1 to 7 do accum:=fpmul(accum, accum, dwords);
		end
		else
			accum:=1;

		if fac>0 then
		begin
			n:=fp2ui(fac);
			fac:=fpipow(fpexp1, n, dwords);
			accum:=accum*fac;
		end;
		if x2<0 then accum:=A0/accum;

		copydf(result, accum, dwords);
	end;

//=======================================================================================


	procedure gauss_leg_rule(const n:Int32; var x: array of decfloat ; var w: array of decfloat; const dwords_in:Int32 = num_dwords);

	VAR
		m, j, i: integer;
		eps, x1, x2, z1, z, xm, xl, pp, p3, p2, p1: decfloat;
		s:string;
	BEGIN
		str(num_digits-5, s);
		s:=trim(s);
		eps:='1e-'+s;
		x1:=-1;
		x2:=1;
		m := (n+1) DIV 2;
		xm := 0.5*(x2+x1);
		xl := 0.5*(x2-x1);
		FOR i := 1 TO m DO BEGIN
			z := cos(3.1415926535897932*(i-0.25)/(n+0.5));
			REPEAT
				p1 := 1;
				p2 := 0;
				FOR j := 1 TO n DO BEGIN
					p3 := p2;
					p2 := p1;
					p1 := ((2*j-1)*z*p2-(j-1)*p3)/j
				END;
				pp := n*(z*p1-p2)/(z*z-1);
				z1 := z;
				z := z1-p1/pp;
			UNTIL (abs(z-z1) <= eps);
			x[i] := xm-xl*z;
			x[n+1-i] := xm+xl*z;
			w[i] := 2*xl/((1-z*z)*pp*pp);
			w[n+1-i] := w[i]
		END
	END;

	function TryStrToDecFloat(Token:String; var Res:decfloat):boolean;
	begin
		res:=str2fp(Token);
		TryStrToDecFloat:=True;
	end;

	function getExponent(n:decfloat): int32;
	var
		ex: int32;
	begin
		If n.exponent<>0 Then
			ex:=(n.exponent And $7FFFFFFF)-BIAS-1
		Else
			ex:=0;
		getExponent:=ex;
	end;

	initialization

	//fppi:=pi_chudnovsky_bs;
	fppi:='3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230665';
	//fpln2:=calcln2;
	fpln2:='0.693147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687542001481020570685734';
	fpPihalf:='1.57079632679489661923132169163975144209858469968755291048747229615390820314310449931401741267105853399107404325664115332';
	fpPi2:='6.28318530717958647692528676655900576839433879875021164194988918461563281257241799725606965068423413596429617302656461329';
	Pi3o2:='4.71238898038468985769396507491925432629575409906265873146241688846172460942931349794205223801317560197322212976992345997';
	fpdtor:='57.2957795130823208767981548141051703324054724665643215491602438612028471483215526324409689958511109441862233816328648933';
	fpPIo6:='0.523598775598298873077107230546583814032861566562517636829157432051302734381034833104672470890352844663691347752213717775';
	fpexp1:='2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992';
	fpSQR3:='1.73205080756887729352744634150587236694280525381038062805580697945193301690880003708114618675724857567562614141540670303';
	fpTMS3:='0.267949192431122706472553658494127633057194746189619371944193020548066983091199962918853813242751424324373858584593296970';
	fpconst1olog10 := '0.434294481903251827651128918916605082294397005803666566114453783165864649208870774729224949338431748318706106744766303734';
	sqrt2pi:='2.50662827463100050241576528481104525300698674060993831662992357634229365460784197494659583837805726611600997266520387964';
	fponec.re :=1; fponec.im := 0;

	gamma_coeff[ 0] := ' 1';
	gamma_coeff[ 1] := ' 0.833333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1';
	gamma_coeff[ 2] := ' 0.347222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222e-2';
	gamma_coeff[ 3] := '-0.268132716049382716049382716049382716049382716049382716049382716049382716049382716049382716049382716049382716049382716049e-2';
	gamma_coeff[ 4] := '-0.229472093621399176954732510288065843621399176954732510288065843621399176954732510288065843621399176954732510288065843621e-3';
	gamma_coeff[ 5] := ' 0.784039221720066627474034881442288849696257103664511071918479325886733294140701548108955516362923770331177738585145992553e-3';
	gamma_coeff[ 6] := ' 0.697281375836585777429398828575783308293596359439980839157793890304178370021991421168375900886188952032573431750386482897e-4';
	gamma_coeff[ 7] := '-0.592166437353693882864836225604401187391585196797816825251667501324565796445083139184648100971832110378063738831914414905e-3';
	gamma_coeff[ 8] := '-0.517179090826059219337057843002058822817853453427298877687537952741427821446111340944445654093573618036352375629924784017e-4';
	gamma_coeff[ 9] := ' 0.839498720672087279993357516764983445198182111593007628503572550042430988763695498779920192645217624472141906814805593786e-3';
	gamma_coeff[10] := ' 0.720489541602001055908571930225015052063451737975478952504489810472170562136200036603978769847471305259792464527192505328e-4';
	gamma_coeff[11] := '-0.191443849856547752650089885832852254487689357895350389611948548607825662367666358433241035717805750334559472625329687275e-2';
	gamma_coeff[12] := '-0.162516262783915816898635123980270998105872591932251887038488056749786172302963502261320355956912097191134273391280343597e-3';
	gamma_coeff[13] := ' 0.640336283380806979482363809026579583040189409396286467630149566540759462780830481022359250933933113060564570086713978430e-2';
	gamma_coeff[14] := ' 0.540164767892604515180467508570241735547254415979217175042122014107572377089417821576193008412964426405197953329177009207e-3';
	gamma_coeff[15] := '-0.295278809456991205054406510546938244465654828254379690562082489809192494566189492309947860518708021787589182390547582114e-1';
	gamma_coeff[16] := '-0.248174360026499773091565836874346432397516804722902067594532201062060654278884665157255284262533712046719408727100189023e-2';
	gamma_coeff[17] := ' .179540117061234856107699407722226330530912823386921731875811924895310683950029446123380319284827094353692100626830387310';
	gamma_coeff[18] := ' 0.150561130400264244123842218771311272602598154554149207382929253316682733257619702959782991646550951436574153935677179582e-1';
	gamma_coeff[19] := '-1.39180109326533748139914776354227314935804561772646246798659253704752404920001483582611648286936847089937525575882720785';
	gamma_coeff[20] := '-.116546276599463200850734036907147969678937334383706618512207226117921453003382408880410622450621616682620424892458577862';
	gamma_coeff[21] := ' 13.3979854551425892176269304320196719504205855650331635741858883489491472699534721421314261069689065942167029506316333619';
	gamma_coeff[22] := ' 1.12080446428991160686263940013992394100874458058978613764384695979903902900510159097655688389529528541466497260708086291';
	gamma_coeff[23] := '-156.801412704022726372823698446041189864295925352586080665153461521199178847908018359693260716382509476562473403701828096';
	gamma_coeff[24] := '-13.1078630226338656590275053222671726562139542672625010472706003948718823093413414154484618522297875238170626571977349659';
	gamma_coeff[25] := ' 2192.55553609052343296901296683540498912174443933789243360552478651203748433247940538907885597303466541735205879449545576';
	gamma_coeff[26] := ' 183.190733484524338088662112060475268304900810166597832590377161208695562878432051985060937330250258211943622604695939859';
	gamma_coeff[27] := '-36101.1192932220759519137910143102123117274408120194227408925427677018309125103754936440090683530295974253049479420862754';
	gamma_coeff[28] := '-3015.07731262230585421582738429513458512616707765571355758280202723181248067955325098556895419104378994426608210212106308';
	gamma_coeff[29] := ' 6.91346376141878121600201494236207859564711767920033294258827498480743934126612121518316826575503168787092449098046130122e5';
	gamma_coeff[30] := ' 57721.3363630407227165872199716323655754083996547324728710542873103200679182342574442942583702988049986042684241190353110';
	gamma_coeff[31] := '-1.52358121457152742068278713804932905941760422599539924300951435704680641286628199952620845994435725967951092409156200895e7';
	gamma_coeff[32] := '-1.27173563925300839222817964702027176992060992077897561600534583165743535363148251322384360656712659249098878060216663536e6';
	gamma_coeff[33] := ' 3.82847679580470263261480955964227610549265280440899139412232695486425990416693311163543854841905047558631260299704454005e8';
	gamma_coeff[34] := ' 3.19498204147977021140035993083213898922678174993726623454135402080918800326622791944558946865120234209162653334685185329e7';
	gamma_coeff[35] := '-1.08809329722618036539009437555168984155254015921379192241333796318914996765708099161421140004353729562174621470822941030e10';
	gamma_coeff[36] := '-9.07894578158421435427828494391453046967517025764055172224013790953916852188602568677468703200571801027349455729606931421e8';
	gamma_coeff[37] := ' 3.47282409152616274523914587858466108775731042219084010629204566179507183135429994192260390419172101208520392823158453465e11';
	gamma_coeff[38] := ' 2.89728458003444433581954284209234687561824957674661343723444683741427313593897887830097155318658141450116989573725366377e10';
	gamma_coeff[39] := '-1.23683936443407959584855563197308673736175627719287502129872485950767770214419271845396015508203187903065994417938345246e13';
	gamma_coeff[40] := '-1.03174023115609124965330104523134253698816723730031128609864915226223664497692184001186566458497814977973891887183579695e12';
	gamma_coeff[41] := ' 4.88745034305486110500979050768819105012419862971381646190382883179981065569248885369401000402792050958590109112728905241e14';
	gamma_coeff[42] := ' 4.07657847965959565639063577711273832105650126536001471723027892792179893974111266769829149244584235877832107770791858656e13';
	gamma_coeff[43] := '-2.13186339168489509268379583288252058181092070724614151850043022458231952550739132100136457830888164022122285454172585535e16';
	gamma_coeff[44] := '-1.77801501138010612594383489135697194994972325216575425915156154102885093037270286462244076491166288792117383150695080832e15';
	gamma_coeff[45] := ' 1.02170115454378600407054121475613221088463288274213614791517561591839912571889361315175429855850647576384766645618922918e18';
	gamma_coeff[46] := ' 8.52054990370652711707338317075144506642892621870230631157040929599356838195883251951682656443684601532311669500126701460e16';
	gamma_coeff[47] := '-5.35719194152387519012295780561959046066966980443756848566863474218620904329642630397343395986538723057423266521787582721e19';
	gamma_coeff[48] := '-4.46737941151948250433339670844108129505566924316203908646197749604896956196494375363913055182427058243715007728564288313e18';
	gamma_coeff[49] := ' 3.06139200177805958311011222548995302823448377890521485292237034053868016292942435146705909580895939245672697755929937441e21';
	gamma_coeff[50] := ' 2.55275989517846398155153280509019670159616634190514288422205265363883284729734355864787738666668462123868550865684162066e20';
	gamma_coeff[51] := '-1.89988531417447383430434877310733370529631362484894834582625350840258649019240280247873481439860829940681036330840040551e23';
	gamma_coeff[52] := '-1.58415162442108182325641963091849383970579359689271719853344369843552150724703215022362935065240658928176022302770473701e22';
	gamma_coeff[53] := ' 1.27627136081367915801215067300665314477682845121638102072079988232338375997771894318275765639190276305818033464552843689e25';
	gamma_coeff[54] := ' 1.06412638579070839964755581940503523265146990642896780639128345161483801313430698890348704227773314532796896717475775501e24';
	gamma_coeff[55] := '-9.25240356525256184896552024400557523282734155133378340830071844736367733804681224799884774089266563770410093584340147328e26';
	gamma_coeff[56] := '-7.71414334093151161077234227538836433275926859518803787051231957971786623043374234077190159351046158671321956308208214381e25';
	gamma_coeff[57] := ' 7.21850102154489005377588591747722010184421130743578095396542188806077207515198294298847714709312468860501041629681583697e28';
	gamma_coeff[58] := ' 6.01817660920955302497359384865412675525992434078231895974986406540249917831074361964465966993740560929553375367380492973e27';
	gamma_coeff[59] := '-6.04493253920648022234407412600966580971414591676908413629403032875627314687846636711347793032394912822256198423228089019e30';
	gamma_coeff[60] := '-5.03959577267195531635704710675864623755563513021623785023406558275111372364348865666981593512745839225915034323626233798e29';
	gamma_coeff[61] := ' 5.42046040306473708910200872754209762441177730163379282757043743442921539110523905079953509795688016001244671642727843737e32';
	gamma_coeff[62] := ' 4.51885202016304859824454532636646620242349131533968002498900697970859511377318343796098082858314174214188083286580131219e31';
	gamma_coeff[63] := '-5.19276945836104974051740979898702938327178803529111899970657164425406241215012798541804204721794107503275200567229286729e34';
	gamma_coeff[64] := '-4.32892308840228496267696927873172028443711667535510556944843702341242872410911147647295403900193969753960042734579149988e33';
	gamma_coeff[65] := ' 5.30347841926986478558769772674374679526277796136877415428148502940282511414845141856441062032111136874747739072762484734e36';
	gamma_coeff[66] := ' 4.42111240092351852910873748610050544762147069511516425082010577171644885907789653023094618489577092720638006104201385432e35';
	gamma_coeff[67] := '-5.76314107379856220649117190736587642225617157185301075869597044116601638378154019530383264664149564585181549421609799558e38';
	gamma_coeff[68] := '-4.80419731202856758575096785773281358475887665368837351388308854722986551856553001836331493820138633663864210224107713758e37';
	gamma_coeff[69] := ' 6.65095547731586233757134602506226681893345991208165369509536504610777944128173027151124148807416763416716096968686517100e40';
	gamma_coeff[70] := ' 5.54417928794743556836580224889720411644919860146941341663154934381022178522660735988852538053152417500498257886159400962e39';

	Bernoulli_fac[ 0] := '-.5';
	Bernoulli_fac[ 1] := '0.083333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333';
	Bernoulli_fac[ 2] := '-0.0013888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888889';
	Bernoulli_fac[ 3] := '3.3068783068783068783068783068783068783068783068783068783068783068783068783068783068783068783068783068783068783068783068783069e-5';
	Bernoulli_fac[ 4] := '-8.2671957671957671957671957671957671957671957671957671957671957671957671957671957671957671957671957671957671957671957671957672e-7';
	Bernoulli_fac[ 5] := '2.0876756987868098979210090321201432312543423654534765645876756987868098979210090321201432312543423654534765645876756987868099e-8';
	Bernoulli_fac[ 6] := '-5.2841901386874931848476822021795566769111742656716201689746663291636836610381583926557471531016504561478106451651425196398741e-10';
	Bernoulli_fac[ 7] := '1.3382536530684678832826980975129123277271425419573567721715869864018012166160314308462456610604758752906901055049203197351345e-11';
	Bernoulli_fac[ 8] := '-3.3896802963225828668301953912494424995721810744115962499612255811031617650561524407359037507393145944862927952460061600013877e-13';
	Bernoulli_fac[ 9] := '8.5860620562778445641359054504256271339539561274865749698784493815783351307329146244656134988702930326189983095847022530973131e-15';
	Bernoulli_fac[10] := '-2.1748686985580618730415164238659178998517916001541660835248118744816299534723419527306624508708524231365576933674021092287274e-16';
	Bernoulli_fac[11] := '5.5090028283602295152026526089022548778615827044688407893777368490753325182933145312127786579586993825947341955861676180549587e-18';
	Bernoulli_fac[12] := '-1.3954464685812523340707686264063549763917636669110995237839659342679375655787601562810342460947433369693884692413053192266794e-19';
	Bernoulli_fac[13] := '3.5347070396294674716932299778037992147245945647149203595069977313091224661427412678677354579054824633933168795522545659677082e-21';
	Bernoulli_fac[14] := '-8.9535174270375468504026113181127410516271392427849625164443408655356959926763835340636411500962239072942372298675938447879001e-23';
	Bernoulli_fac[15] := '2.2679524523376830603109507388681660632203543297441712610647033268678298964148063678381848536080108864566115725797124055187523e-24';
	Bernoulli_fac[16] := '-5.7447906688722024452638819876070183996247766917890951717645030435378005100088496029056318139535163521779370185248190239187365e-26';
	Bernoulli_fac[17] := '1.4551724756148649018662648672713293357208889558290473512886334021214778052178986025990907828880705045801862118731872549000038e-27';
	Bernoulli_fac[18] := '-3.6859949406653101781817824799086603744462982064359608441829046575546291278460091805256186629695399222974516768303680768919263e-29';
	Bernoulli_fac[19] := '9.3367342570950446720325551527856232954436887120383714258605088611901188869432090607490011453277511709410363960088353129119055e-31';
	Bernoulli_fac[20] := '-2.3650224157006299345596351963698382400696562503497727951255301437351689313234761154145295718127071061761411293926243451325562e-32';
	Bernoulli_fac[21] := '5.9906717624821343046599123968196578264493699903956181329666424130797629471793519981612323839614554852568239040806726588995128e-34';
	Bernoulli_fac[22] := '-1.5174548844682902617108131358647189315408843012454426869736145618143953170448507246056796034899754184023228866918414497661058e-35';
	Bernoulli_fac[23] := '3.8437581254541882322294452909902321059018090473728865398130796252375142545635062229304064951688992943554191150034177620955737e-37';
	Bernoulli_fac[24] := '-9.7363530726466910352676212792504541809551090795412416887734353583087034631549823547563402450384308326078748249755023205494004e-39';
	Bernoulli_fac[25] := '2.4662470442006809571064002802888428859241773384046221997933670889265767289944168994090722427737443266392054145488517912442750e-40';
	Bernoulli_fac[26] := '-6.2470767418207436931487567947233686925765739802212621514703689323745743093235073453685443080489320400409172441699517721376355e-42';
	Bernoulli_fac[27] := '1.5824030244644914297510817068287639403286027624179166162140253039015117418999973525802084212601446642620014330093351581482364e-43';
	Bernoulli_fac[28] := '-4.0082736859489359685300121905219826626811254808496040277438964116592033491787091000561842595652054045204200928532874915711740e-45';
	Bernoulli_fac[29] := '1.0153075855569556311630713945378762327067794638832996356827823841161264078193143189566751697785126763094589814321700842332273e-46';
	Bernoulli_fac[30] := '-2.5718041582418717499248194097644548855731577502217601608187158823068058926052257241623481922798960208725987982293323567979923e-48';
	Bernoulli_fac[31] := '6.5144560352338149315584348586418580231420396480615062396978476728376073823296528061321665580121493992109263264830561135065583e-50';
	Bernoulli_fac[32] := '-1.6501309906896524555060987804793230091879183091243903605141108515641747255867244232697460463290438370921231788315688965726184e-51';
	Bernoulli_fac[33] := '4.1798306285394758948501872347094070329312870995339112249216314280651132794146888010437985296100408748943807374825010631560910e-53';
	Bernoulli_fac[34] := '-1.0587634667702908770270420242791172873352071246178360650733876162286683365066279055055624176846828793628136640603977895140372e-54';
	Bernoulli_fac[35] := '2.6818791912607706661409848588415103397690942686090118518199045742904944496379348032175308234124548440609863557261262067203289e-56';
	Bernoulli_fac[36] := '-6.7932793511074212095271802995338946118945424532313411277912191561000857066486792271304987266790652382941286406536472063521086e-58';
	Bernoulli_fac[37] := '1.7207577616681404905363499407582306642815875179953978007837551911155802976989232688268053450129336488306759824512671274944218e-59';
	Bernoulli_fac[38] := '-4.3587303293488938434001998497731611091275265226574759731848503961372269004550512418303885289132772666097081965529139565453541e-61';
	Bernoulli_fac[39] := '1.1040792903684666750838395976444273230873858069236861903982817911816599858283442448029326267690780890267138648327870514700816e-62';
	Bernoulli_fac[40] := '-2.7966655133781345072047937531186265538644898338925587556702061386797353780378035840064465002672448127140392805384426001557633e-64';
	Bernoulli_fac[41] := '7.0840365016794701985093884223803493343233878282522421642099970759945234842557107413733957215035945077422862362604633134885951e-66';
	Bernoulli_fac[42] := '-1.7944074082892240666052573093367559867049219739052200233591338702660539989433752258121897666477668900423608741864577248701312e-67';
	Bernoulli_fac[43] := '4.5452870636110961070850791246394481858634198268973747159896656439073580551312677350016771929489373538799660061307800482311061e-69';
	Bernoulli_fac[44] := '-1.1513346631982051812730029007862433774602680136860452024051019678942818548598882354218764600351190069679690497976404891059352e-70';
	Bernoulli_fac[45] := '2.9163647710923613547033689800529707427867341575479118965689807463964968836710042873107629944417727479611388316499039003266602e-72';
	Bernoulli_fac[46] := '-7.3872382634973375625733753947097366977865048567538833277733827905496569232612305170338969369969740825883307200675668667906483e-74';
	Bernoulli_fac[47] := '1.8712093117637953062252618714080746593010188462223387105749305312506675209198976490468436789802816315492716864340430793564961e-75';
	Bernoulli_fac[48] := '-4.7398285577617994054995634412108870791977492250801906167932740135781420207257486562916149301494430186154883378407327817895050e-77';
	Bernoulli_fac[49] := '1.2006125993354506519817100221141291236516247116986875074894710471333276473097521905209014196034240306126510157356033760058130e-78';
	Bernoulli_fac[50] := '-3.0411872415142923830371220748388955327985525575352100010894014134687189922768475534473571677972659869760016329063028357636522e-80';
	Bernoulli_fac[51] := '7.7034172747051062728765033945972188824524480656164977385223628350290194350977106782537710081174931940019763876882236045056084e-82';
	Bernoulli_fac[52] := '-1.9512983909098830711123237681456177474867779326287312595522051132966077288751382333729832893573179473085136977158156347658726e-83';
	Bernoulli_fac[53] := '4.9426965651594614748964001538883279702750864876010973614220209218549491715347721254941082092007349235608435885433357725849195e-85';
	Bernoulli_fac[54] := '-1.2519996659171847921681857688255382256613101478412675018592588178828836457191166902122390118148211406553818926277126960219954e-86';
	Bernoulli_fac[55] := '3.1713522017635154606453951675914196545944083062098027978712715309801652602710635756367678047394069537203930193795779584514166e-88';
	Bernoulli_fac[56] := '-8.0331289707353344613949571241406043209088638963578510073763936443200008085813970439838379786607691377351355750344671512984656e-90';
	Bernoulli_fac[57] := '2.0348153391661465707815120758715955705475032423806398081830833625138414563757764533056587450495583972919433964460588577701957e-91';
	Bernoulli_fac[58] := '-5.1542474664474738590695091367793757748203475533416550414602754469758875009672963116699001406091368840898943186791750701879607e-93';
	Bernoulli_fac[59] := '1.3055861352149467245765754771001878937795610150492203169516353232547309748059062099844996547108863817766847011194764937493260e-94';
	Bernoulli_fac[60] := '-3.3070883141750912484517539812530574149704848311262370496984442173454728473149791235560687984648602964061482523358446895707170e-96';
	Bernoulli_fac[61] := '8.3769525600490913030350758435207266985008710064992915741986034103423156429232264354728364364404641457474303932629014493099994e-98';
	Bernoulli_fac[62] := '-2.1219068717497137695288988470380586224930360752040247151031858267734232587429554172609862213504937658163433086516958175673080e-99';
	Bernoulli_fac[63] := '5.3748528956122802562805569922432021936713624712660056831762697398539258073946333524686784703862941350720434938271905157373473e-101';
	Bernoulli_fac[64] := '-1.3614661432172069392504095313379850565511866355903681931540330732801952476146447849527554330335635628811653672697376944778833e-102';
	Bernoulli_fac[65] := '3.4486340279933990342773740495738169045523462689026039422925143149193185157658126589134674750139047308533196997802530467593047e-104';
	Bernoulli_fac[66] := '-8.7354920416383550602712612095027905086745403529670836331150970649504379863583792027487786805512931521742833099477219344479890e-106';
	Bernoulli_fac[67] := '2.2127259833925497084311154562133174043197738400769687003309471582495933269591405957888650416139821756047196774128163492895982e-107';
	Bernoulli_fac[68] := '-5.6049003928372241708220677872990576941581603552833772169514808702876377589775618950228520161628852016026330094604637844383114e-109';
	Bernoulli_fac[69] := '1.4197378549991787673096827703627044937522866246369183559498810588789733295140708701154925117312331819617562374464955930216789e-110';
	Bernoulli_fac[70] := '-3.5962379982587626636746239711915365171205632677582896074576877036203126876141204842243660789970562591295585940248627969593212e-112';
	Bernoulli_fac[71] := '9.1093772660782318645768927671397510232496074791943252826363167388985839592405108851705809034137670417908698807050322736385344e-114';
	Bernoulli_fac[72] := '-2.3074322171091232885036435704648503002329764148062133408862694692557108612514330596154043037255634752776031335467914420632543e-115';
	Bernoulli_fac[73] := '5.8447940852990019944930889298673049324080320862745050892736089565998959017253113877436706809946961577787694427723173071665301e-117';
	Bernoulli_fac[74] := '-1.4805036371705744952518203162987911250812478035503350461762093493242144958402456005809476943846261259612242824335728146559890e-118';
	Bernoulli_fac[75] := '3.7501595226227196890811761291914005623496889627156492621381709846159226711084583847482690230683105319625039092180705894514217e-120';
	
end.
