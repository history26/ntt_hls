/*
 * Copyright (C) 2022 DZK
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _ARITHMETIC_FUNC_H_
#define _ARITHMETIC_FUNC_H_

#include <ap_int.h>

//R=b^n,b=2^LONG_WIDTH,n=1
#define BASE_WIDTH 256
#define in_WIDTH 512
typedef ap_uint<in_WIDTH> in_type;
typedef ap_uint<257> long_uint;
typedef ap_uint<512> long_long_uint;
typedef ap_int<512> long_int;
typedef ap_uint<32> IDXTYPE;

//p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
const ap_uint<257> p("0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001",16);

/// INV = -(p^{-1} mod 2^64) mod 2^6
// const long_uint INV= 0x992d30ecffffffff;
// base = 2^256
// const long_uint INV = 0x1E0E3B00E1DD872A 8F34D6691037659A 75A6DE91C8D4FCC3 C9EDA265AC589658;
const ap_uint<257> INV("0x1E0E3B00E1DD872A8F34D6691037659A75A6DE91C8D4FCC3C9EDA265AC589658",16);

long_uint add_long(long_uint left,long_uint right);
long_uint sub_long(long_uint left,long_uint right);
long_long_uint mult_long(long_long_uint left,long_long_uint right);
long_uint mult_mod_base_long(long_long_uint left,long_long_uint right);
long_uint mul_partly(long_uint left,long_uint right);
long_long_uint mul_long_opt(long_uint left,long_uint right);
long_uint to_inverse(long_uint p,long_uint R);
long_uint montgomery_reduce(long_long_uint t,long_uint m,long_uint negat_m_inv);
//res = left * right ,then montgomery_reduce
long_uint MUL_MOD(long_long_uint left,long_long_uint right,long_uint m);
long_uint SUB_MOD(long_uint left,long_uint right,long_uint m);
long_uint ADD_MOD(long_uint left,long_uint right,long_uint m);

const long_uint EPSILON = long_uint((1ull << 32) - 1);
const long_uint ZERO = long_uint(0);
const long_uint ONE = long_uint(1);
const long_uint TWO = long_uint(2);

//18446744069414584321
//const GF MODULUS = GF(0xFFFFFFFF00000001);

const long_uint OMEGA[25]={
		ap_uint<257>("0x0001",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<257>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16)
	};

//const GF OMEGA[33] = {
//	GF(1ull),						// for a domain of 2^0
//	GF(18446744069414584320ull),	// for a domain of 2^1
//	GF(281474976710656ull),			// for a domain of 2^2
//	GF(18446744069397807105ull),	// for a domain of 2^3
//	GF(17293822564807737345ull),	// for a domain of 2^4
//	GF(70368744161280ull),			// for a domain of 2^5
//	GF(549755813888ull),			// for a domain of 2^6
//	GF(17870292113338400769ull),	// for a domain of 2^7
//	GF(13797081185216407910ull),	// for a domain of 2^8
//	GF(1803076106186727246ull),		// for a domain of 2^9
//	GF(11353340290879379826ull),	// for a domain of 2^10
//	GF(455906449640507599ull),		// for a domain of 2^11
//	GF(17492915097719143606ull),	// for a domain of 2^12
//	GF(1532612707718625687ull),		// for a domain of 2^13
//	GF(16207902636198568418ull),	// for a domain of 2^14
//	GF(17776499369601055404ull),	// for a domain of 2^15
//	GF(6115771955107415310ull),		// for a domain of 2^16
//	GF(12380578893860276750ull),	// for a domain of 2^17
//	GF(9306717745644682924ull),		// for a domain of 2^18
//	GF(18146160046829613826ull),	// for a domain of 2^19
//	GF(3511170319078647661ull),		// for a domain of 2^20
//	GF(17654865857378133588ull),	// for a domain of 2^21
//	GF(5416168637041100469ull),		// for a domain of 2^22
//	GF(16905767614792059275ull),	// for a domain of 2^23
//	GF(9713644485405565297ull),		// for a domain of 2^24
//	GF(5456943929260765144ull),		// for a domain of 2^25
//	GF(17096174751763063430ull),	// for a domain of 2^26
//	GF(1213594585890690845ull), 	// for a domain of 2^27
//	GF(6414415596519834757ull), 	// for a domain of 2^28
//	GF(16116352524544190054ull), 	// for a domain of 2^29
//	GF(9123114210336311365ull), 	// for a domain of 2^30
//	GF(4614640910117430873ull), 	// for a domain of 2^31
//	GF(1753635133440165772ull) 		// for a domain of 2^32
//};


#endif
