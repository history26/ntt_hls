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

#ifndef _NTT_FUNC_H_
#define _NTT_FUNC_H_

#include "arithmetic.hpp"
#include "define.hpp"
#include <iostream>


//void ntt_core(GF w, GF operanda_in, GF operandb_in, GF &operanda_out, GF &operandb_out);
//void ntt_2_12_4core(GF in[bramnum][bramsize]);

void ntt_core_large_bw(long_uint w,long_uint m, long_uint operand_ina,long_uint operand_inb,
						long_uint &operand_outa,long_uint &operand_outb);
void NTT_2_12_in_place_large_bw(long_uint in[4096]);
void ntt_2_12_2core_large_bw(long_uint in[bramnum][bramsize]);
void ntt_2_12_4core_large_bw(long_uint in[bramnum][bramsize]);
void ntt_2_12_8core_large_bw(long_uint in[bramnum][bramsize]);
void ntt_2_12_16core_large_bw(long_uint in[bramnum][bramsize]);


#endif
