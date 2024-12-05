#ifndef _NTT_V2_FUNC_H_
#define _NTT_V2_FUNC_H_

#include "arithmetic.hpp"
#include "define.hpp"
#include <iostream>
#include "ntt.hpp"

void ntt_4core_new(long_uint operand[bramnum][bramsize],
				   long_uint rp[RPBRAMNUM][RPBRAMSIZE],
				   long_uint modulus);
void ntt_4core_new_v2(long_uint operand[bramnum][bramsize],
				   long_uint rp[RPBRAMNUM][RPBRAMSIZE],
				   long_uint modulus);

#endif
