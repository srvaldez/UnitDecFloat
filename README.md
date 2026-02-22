# UnitDecFloat
decimal floating-point unit for FreePascal


the unit overloads the arithmetic and logical operators as well as the common mathematical functions including the zeta and gamma functions, all functions support the real and complex type


the unit uses the following records

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
each dword in the mantissa holds 8 digits


the math functions are tailored for about 100+ digits of precision, however the basic arithmetic functions and the square root function are only limited by the precision of the mantissa
