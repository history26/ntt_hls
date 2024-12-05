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

#include "ntt.hpp"
using namespace std;


void ntt_core_large_bw(long_uint w,long_uint m, long_uint operand_ina,long_uint operand_inb,
						long_uint &operand_outa,long_uint &operand_outb){
	#pragma HLS INLINE
	operand_inb = MUL_MOD(operand_inb, w, m);
	operand_outa = ADD_MOD(operand_ina, operand_inb, m);
	operand_outb = SUB_MOD(operand_ina, operand_inb, m);
}



void NTT_2_12_in_place_large_bw(long_uint in[4096]) {
	const int n = 4096;
	const int logn = 12;

	typedef ap_uint<logn> INDEX;

	long_uint tmp;
	for (INDEX k = INDEX(0);;) {
		INDEX rk = k;
		rk.reverse();
		if (k < rk) {
			tmp = in[rk];
			in[rk] = in[k];
			in[k] = tmp;
		}

		k++;
		if (k.iszero()) { // overflow
			break;
		}
	}

	// // for(int i=0;i<4096;i++){
	// // 	if(in[i]==2){
	// // 		cout<<"i = "<<i<<" data in is "<<in[i]<<endl;
	// // 	}

	// }

	int m = 1;
	for (int i = 1; i <= logn; i++) {
		// w_m is 2^i-th root of unity
		long_uint w_m = OMEGA[i];

		int k = 0;
		while (k < n) {
			// w = w_m^j at the start of every loop iteration
			long_uint w = 1;

			for (int j = 0; j < m; j++) {
				long_uint t = in[k + j + m];
				t = MUL_MOD(t, w, p);

				long_uint tmp = in[k + j];
				tmp = SUB_MOD(tmp, t, p);

				in[k + j + m] = tmp;
				in[k + j] = ADD_MOD(in[k + j], t, p);

				w = MUL_MOD(w, w_m, p);
			}

			k += 2 * m;
		}

		m *= 2;
	}
}


void ntt_2_12_2core_large_bw(long_uint in[bramnum][bramsize]) {
#pragma HLS BIND_STORAGE variable=OMEGA type=rom_np impl=lutram

	const int n = 4096;
	const int logn = 12;

	uint32_t stepsize;
	uint32_t stepnum;
	uint32_t MEa, MEb;
	long_uint w[corenum], w_m;
	long_uint operanda[corenum], operandb[corenum], operandaout[corenum], operandbout[corenum];
	long_uint operanda_[corenum], operandb_[corenum], operandaout_[corenum], operandbout_[corenum];

	ntt_stage1:for(uint32_t i = 1; i <= L_BRAMSIZE; i++) {
		w[0] = 1; // no need to use w[1] in ntt_stage1
		w_m = OMEGA[i]; // w_m is 2^i-th root of unity
		stepsize = 1 << (i-1);
		stepnum = bramsize;
		MEa = 0;
		MEb = stepsize;

		ntt_stage1_inner:for(uint32_t j = 0; j < stepnum; j+=2){
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = in inter false
			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				operanda[m] = in[m][MEa];
				operandb[m] = in[m][MEb];
				operanda_[m] = in[m + corenum][MEa];
				operandb_[m] = in[m + corenum][MEb];
			}

			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				// ntt core
				ntt_core_large_bw(w[0], p, operanda[m], operandb[m], operandaout[m], operandbout[m]);
				ntt_core_large_bw(w[0], p, operanda_[m], operandb_[m], operandaout_[m], operandbout_[m]);
			}

			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				in[m][MEa] = operandaout[m];
				in[m][MEb] = operandbout[m];
				in[m + corenum][MEa] = operandaout_[m];
				in[m + corenum][MEb] = operandbout_[m];
			}

			if(((MEa + 1) & stepsize)){
				MEa += stepsize + 1;
				MEb += stepsize + 1;
				w[0] = ONE;
			}
			else{
				MEa += 1;
				MEb += 1;
				w[0] = MUL_MOD(w[0], w_m, p);  // compute w for next ntt_core
			}
		}
	}

	ntt_stage2:for(uint32_t i = L_BRAMSIZE + 1; i <= logn; i++) {
			w[0] = 1;
			w[1] = (i==12)?OMEGA[2]:ONE;  // the w value is precomputed, which is the value of w when j=1024
			w_m = OMEGA[i]; // w_m is 2^i-th root of unity
			stepsize = 1 << (i-1);
			stepnum = bramsize;
			MEa = 0;
			MEb = 0;

			ntt_stage2_inner:for(uint32_t j = 0; j < stepnum; j++){
#pragma HLS PIPELINE II = 1
#pragma HLS DEPENDENCE variable = in inter false
				switch(i){
					case 11:
						operanda[0] = in[0][MEa];
						operandb[0] = in[1][MEb];
						operanda[1] = in[2][MEa];
						operandb[1] = in[3][MEb];
						break;
					case 12:
						operanda[0] = in[0][MEa];
						operandb[0] = in[2][MEb];
						operanda[1] = in[1][MEa];
						operandb[1] = in[3][MEb];
				}

				for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
					// ntt core
					ntt_core_large_bw(w[m], p, operanda[m], operandb[m], operandaout[m], operandbout[m]);
				}


				switch(i){
					case 11:
						in[0][MEa] = operandaout[0];
						in[1][MEb] = operandbout[0];
						in[2][MEa] = operandaout[1];
						in[3][MEb] = operandbout[1];
						break;
					case 12:
						in[0][MEa] = operandaout[0];
						in[2][MEb] = operandbout[0];
						in[1][MEa] = operandaout[1];
						in[3][MEb] = operandbout[1];
				}


				MEa += 1;
				MEb += 1;

				for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
					w[m] = MUL_MOD(w[m], w_m, p);  // compute w for next ntt_core
				}


			}
		}
}

void ntt_2_12_4core_large_bw(long_uint in[bramnum][bramsize]) {
#pragma HLS BIND_STORAGE variable=OMEGA type=rom_np impl=lutram

	const int n = 4096;
	const int logn = 12;

	uint32_t stepsize;
	uint32_t stepnum;
	int MEa, MEb;
	long_uint w[corenum], w_m;
	long_uint operanda[corenum], operandb[corenum], operandaout[corenum], operandbout[corenum];
	long_uint operanda_[corenum], operandb_[corenum], operandaout_[corenum], operandbout_[corenum];

	ntt_stage1:for(uint32_t i = 1; i <= L_BRAMSIZE; i++) {
		w[0] = ONE; // no need to use w[1], in ntt_stage1
		w_m = OMEGA[i]; // w_m is 2^i-th root of unity
		stepsize = 1 << (i-1);
		stepnum = bramsize;
		MEa = 0;
		MEb = stepsize;

		ntt_stage1_inner:for(uint32_t j = 0; j < stepnum; j+=2){
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = in inter false
			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				operanda[m] = in[m][MEa];
				operandb[m] = in[m][MEb];
				operanda_[m] = in[m + corenum][MEa];
				operandb_[m] = in[m + corenum][MEb];
			}

			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				// ntt core
				ntt_core_large_bw(w[0],p, operanda[m], operandb[m], operandaout[m], operandbout[m]);
				ntt_core_large_bw(w[0],p, operanda_[m], operandb_[m], operandaout_[m], operandbout_[m]);
			}

			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				in[m][MEa] = operandaout[m];
				in[m][MEb] = operandbout[m];
				in[m + corenum][MEa] = operandaout_[m];
				in[m + corenum][MEb] = operandbout_[m];
			}

			if(((MEa + 1) & stepsize)){
				MEa += stepsize + 1;
				MEb += stepsize + 1;
				w[0] = ONE;
			}
			else{
				MEa += 1;
				MEb += 1;
				w[0] = MUL_MOD(w[0], w_m,p);  // compute w for next ntt_core
			}
		}
	}

	ntt_stage2:for(uint32_t i = L_BRAMSIZE + 1; i <= logn; i++) {
			w[0] = ONE;
			w_m = OMEGA[i]; // w_m is 2^i-th root of unity
			stepsize = 1 << (i-1);
			stepnum = bramsize;
			MEa = 0;
			MEb = 0;
			ntt_stage2_inner:for(uint32_t j = 0; j < stepnum; j++){
#pragma HLS PIPELINE II = 1
#pragma HLS DEPENDENCE variable = in inter false
				switch(i){
					case 10:
						operanda[0] = in[0][MEa];
						operandb[0] = in[1][MEb];
						operanda[1] = in[2][MEa];
						operandb[1] = in[3][MEb];
						operanda[2] = in[4][MEa];
						operandb[2] = in[5][MEb];
						operanda[3] = in[6][MEa];
						operandb[3] = in[7][MEb];
						w[1] = ONE;
						w[2] = ONE;
						w[3] = ONE;
						break;
					case 11:
						operanda[0] = in[0][MEa];
						operandb[0] = in[2][MEb];
						operanda[1] = in[1][MEa];
						operandb[1] = in[3][MEb];
						operanda[2] = in[4][MEa];
						operandb[2] = in[6][MEb];
						operanda[3] = in[5][MEa];
						operandb[3] = in[7][MEb];
						w[1] = OMEGA[2];
						w[2] = ONE;
						w[3] = OMEGA[2];
						break;
					case 12:
						operanda[0] = in[0][MEa];
						operandb[0] = in[4][MEb];
						operanda[1] = in[1][MEa];
						operandb[1] = in[5][MEb];
						operanda[2] = in[2][MEa];
						operandb[2] = in[6][MEb];
						operanda[3] = in[3][MEa];
						operandb[3] = in[7][MEb];
						w[1] = OMEGA[3];
						w[2] = OMEGA[2];
						w[3] = long_uint("18446742969902956801",10);
				}

				for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
					// ntt core
					ntt_core_large_bw(w[m],p, operanda[m], operandb[m], operandaout[m], operandbout[m]);
				}


				switch(i){
					case 10:
						in[0][MEa] = operandaout[0];
						in[1][MEb] = operandbout[0];
						in[2][MEa] = operandaout[1];
						in[3][MEb] = operandbout[1];
						in[4][MEa] = operandaout[2];
						in[5][MEb] = operandbout[2];
						in[6][MEa] = operandaout[3];
						in[7][MEb] = operandbout[3];
						break;
					case 11:
						in[0][MEa] = operandaout[0];
						in[2][MEb] = operandbout[0];
						in[1][MEa] = operandaout[1];
						in[3][MEb] = operandbout[1];
						in[4][MEa] = operandaout[2];
						in[6][MEb] = operandbout[2];
						in[5][MEa] = operandaout[3];
						in[7][MEb] = operandbout[3];
						break;
					case 12:
						in[0][MEa] = operandaout[0];
						in[4][MEb] = operandbout[0];
						in[1][MEa] = operandaout[1];
						in[5][MEb] = operandbout[1];
						in[2][MEa] = operandaout[2];
						in[6][MEb] = operandbout[2];
						in[3][MEa] = operandaout[3];
						in[7][MEb] = operandbout[3];
				}


				MEa += 1;
				MEb += 1;

				for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
					w[m] = MUL_MOD(w[m], w_m,p);  // compute w for next ntt_core
				}


			}
		}
}


void ntt_2_12_8core_large_bw(long_uint in[bramnum][bramsize]) {
#pragma HLS BIND_STORAGE variable=OMEGA type=rom_np impl=lutram

	const int n = 4096;
	const int logn = 12;

	uint32_t stepsize;
	uint32_t stepnum;
	uint32_t MEa, MEb;
	long_uint w[corenum], w_m;
	long_uint operanda[corenum], operandb[corenum], operandaout[corenum], operandbout[corenum];
	long_uint operanda_[corenum], operandb_[corenum], operandaout_[corenum], operandbout_[corenum];

	ntt_stage1:for(uint32_t i = 1; i <= L_BRAMSIZE; i++) {
		w[0] = ONE; // no need to use w[1], in ntt_stage1
		w_m = OMEGA[i]; // w_m is 2^i-th root of unity
		stepsize = 1 << (i-1);
		stepnum = bramsize;
		MEa = 0;
		MEb = stepsize;

		ntt_stage1_inner:for(uint32_t j = 0; j < stepnum; j+=2){
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = in inter false
			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				operanda[m] = in[m][MEa];
				operandb[m] = in[m][MEb];
				operanda_[m] = in[m + corenum][MEa];
				operandb_[m] = in[m + corenum][MEb];
			}

			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				// ntt core
				ntt_core_large_bw(w[0],p, operanda[m], operandb[m], operandaout[m], operandbout[m]);
				ntt_core_large_bw(w[0],p, operanda_[m], operandb_[m], operandaout_[m], operandbout_[m]);
			}

			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				in[m][MEa] = operandaout[m];
				in[m][MEb] = operandbout[m];
				in[m + corenum][MEa] = operandaout_[m];
				in[m + corenum][MEb] = operandbout_[m];
			}

			if(((MEa + 1) & stepsize)){
				MEa += stepsize + 1;
				MEb += stepsize + 1;
				w[0] = ONE;
			}
			else{
				MEa += 1;
				MEb += 1;
				w[0] = MUL_MOD(w[0], w_m,p);  // compute w for next ntt_core
			}
		}
	}

	ntt_stage2:for(uint32_t i = L_BRAMSIZE + 1; i <= logn; i++) {
		//precompute w
		switch(i){
			case 9:
				w[0] = ONE;
				w[1] = ONE;
				w[2] = ONE;
				w[3] = ONE;
				w[4] = ONE;
				w[5] = ONE;
				w[6] = ONE;
				w[7] = ONE;
				break;
			case 10:
				w[0] = ONE;
				w[1] = OMEGA[2];
				w[2] = ONE;
				w[3] = OMEGA[2];
				w[4] = ONE;
				w[5] = OMEGA[2];
				w[6] = ONE;
				w[7] = OMEGA[2];
				break;
			case 11:
				w[0] = ONE;
				w[1] = OMEGA[3];
				w[1] = OMEGA[3];
				w[2] = OMEGA[2];
				w[3] = long_uint("18446742969902956801",10);
				w[4] = ONE;
				w[5] = OMEGA[3];
				w[6] = OMEGA[2];
				w[7] = long_uint("18446742969902956801",10);
				break;
			case 12:
				w[0] = ONE;
				w[1] = OMEGA[4]; // j=256
				w[2] = OMEGA[3]; // j=512
				w[3] = long_uint("4503599626321920",10);     // j=768
				w[4] = OMEGA[2];      // j=1024
				w[5] = long_uint("18446744069414588417",10); // j=1280
				w[6] = long_uint("18446742969902956801",10); // j=1536
				w[7] = long_uint("18446744000695107585",10); // j=1792
		}
		w_m = OMEGA[i]; // w_m is 2^i-th root of unity
		stepsize = 1 << (i-1);
		stepnum = bramsize;
		MEa = 0;
		MEb = 0;
		ntt_stage2_inner:for(uint32_t j = 0; j < stepnum; j++){
#pragma HLS PIPELINE II = 1
#pragma HLS DEPENDENCE variable = in inter false
			switch(i){
				case 9:
					operanda[0] = in[0][MEa];
					operandb[0] = in[1][MEb];
					operanda[1] = in[2][MEa];
					operandb[1] = in[3][MEb];
					operanda[2] = in[4][MEa];
					operandb[2] = in[5][MEb];
					operanda[3] = in[6][MEa];
					operandb[3] = in[7][MEb];
					operanda[4] = in[8][MEa];
					operandb[4] = in[9][MEb];
					operanda[5] = in[10][MEa];
					operandb[5] = in[11][MEb];
					operanda[6] = in[12][MEa];
					operandb[6] = in[13][MEb];
					operanda[7] = in[14][MEa];
					operandb[7] = in[15][MEb];
					break;
				case 10:
					operanda[0] = in[0][MEa];
					operandb[0] = in[2][MEb];
					operanda[1] = in[1][MEa];
					operandb[1] = in[3][MEb];
					operanda[2] = in[4][MEa];
					operandb[2] = in[6][MEb];
					operanda[3] = in[5][MEa];
					operandb[3] = in[7][MEb];
					operanda[4] = in[8][MEa];
					operandb[4] = in[10][MEb];
					operanda[5] = in[9][MEa];
					operandb[5] = in[11][MEb];
					operanda[6] = in[12][MEa];
					operandb[6] = in[14][MEb];
					operanda[7] = in[13][MEa];
					operandb[7] = in[15][MEb];
					break;
				case 11:
					operanda[0] = in[0][MEa];
					operandb[0] = in[4][MEb];
					operanda[1] = in[1][MEa];
					operandb[1] = in[5][MEb];
					operanda[2] = in[2][MEa];
					operandb[2] = in[6][MEb];
					operanda[3] = in[3][MEa];
					operandb[3] = in[7][MEb];
					operanda[4] = in[8][MEa];
					operandb[4] = in[12][MEb];
					operanda[5] = in[9][MEa];
					operandb[5] = in[13][MEb];
					operanda[6] = in[10][MEa];
					operandb[6] = in[14][MEb];
					operanda[7] = in[11][MEa];
					operandb[7] = in[15][MEb];
					break;
				case 12:
					operanda[0] = in[0][MEa];
					operandb[0] = in[8][MEb];
					operanda[1] = in[1][MEa];
					operandb[1] = in[9][MEb];
					operanda[2] = in[2][MEa];
					operandb[2] = in[10][MEb];
					operanda[3] = in[3][MEa];
					operandb[3] = in[11][MEb];
					operanda[4] = in[4][MEa];
					operandb[4] = in[12][MEb];
					operanda[5] = in[5][MEa];
					operandb[5] = in[13][MEb];
					operanda[6] = in[6][MEa];
					operandb[6] = in[14][MEb];
					operanda[7] = in[7][MEa];
					operandb[7] = in[15][MEb];
			}

			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				// ntt core
				ntt_core_large_bw(w[m],p, operanda[m], operandb[m], operandaout[m], operandbout[m]);
			}


			switch(i){
				case 9:
					in[0][MEa] = operandaout[0];
					in[1][MEb] = operandbout[0];
					in[2][MEa] = operandaout[1];
					in[3][MEb] = operandbout[1];
					in[4][MEa] = operandaout[2];
					in[5][MEb] = operandbout[2];
					in[6][MEa] = operandaout[3];
					in[7][MEb] = operandbout[3];
					in[8][MEa] = operandaout[4];
					in[9][MEb] = operandbout[4];
					in[10][MEa] = operandaout[5];
					in[11][MEb] = operandbout[5];
					in[12][MEa] = operandaout[6];
					in[13][MEb] = operandbout[6];
					in[14][MEa] = operandaout[7];
					in[15][MEb] = operandbout[7];
					break;
				case 10:
					in[0][MEa] = operandaout[0];
					in[2][MEb] = operandbout[0];
					in[1][MEa] = operandaout[1];
					in[3][MEb] = operandbout[1];
					in[4][MEa] = operandaout[2];
					in[6][MEb] = operandbout[2];
					in[5][MEa] = operandaout[3];
					in[7][MEb] = operandbout[3];
					in[8][MEa] = operandaout[4];
					in[10][MEb] = operandbout[4];
					in[9][MEa] = operandaout[5];
					in[11][MEb] = operandbout[5];
					in[12][MEa] = operandaout[6];
					in[14][MEb] = operandbout[6];
					in[13][MEa] = operandaout[7];
					in[15][MEb] = operandbout[7];
					break;
				case 11:
					in[0][MEa] = operandaout[0];
					in[4][MEb] = operandbout[0];
					in[1][MEa] = operandaout[1];
					in[5][MEb] = operandbout[1];
					in[2][MEa] = operandaout[2];
					in[6][MEb] = operandbout[2];
					in[3][MEa] = operandaout[3];
					in[7][MEb] = operandbout[3];
					in[8][MEa] = operandaout[4];
					in[12][MEb] = operandbout[4];
					in[9][MEa] = operandaout[5];
					in[13][MEb] = operandbout[5];
					in[10][MEa] = operandaout[6];
					in[14][MEb] = operandbout[6];
					in[11][MEa] = operandaout[7];
					in[15][MEb] = operandbout[7];
					break;
				case 12:
					in[0][MEa] = operandaout[0];
					in[8][MEb] = operandbout[0];
					in[1][MEa] = operandaout[1];
					in[9][MEb] = operandbout[1];
					in[2][MEa] = operandaout[2];
					in[10][MEb] = operandbout[2];
					in[3][MEa] = operandaout[3];
					in[11][MEb] = operandbout[3];
					in[4][MEa] = operandaout[4];
					in[12][MEb] = operandbout[4];
					in[5][MEa] = operandaout[5];
					in[13][MEb] = operandbout[5];
					in[6][MEa] = operandaout[6];
					in[14][MEb] = operandbout[6];
					in[7][MEa] = operandaout[7];
					in[15][MEb] = operandbout[7];
			}


			MEa += 1;
			MEb += 1;

			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				w[m] = MUL_MOD(w[m], w_m,p);  // compute w for next ntt_core
			}


		}
	}
}



void ntt_2_12_16core_large_bw(long_uint in[bramnum][bramsize])
{

#pragma HLS BIND_STORAGE variable=OMEGA type=rom_np impl=lutram

	const int n = 4096;
	const int logn = 12;

	uint32_t stepsize;
	uint32_t stepnum;
	uint32_t MEa, MEb;
	long_uint w[corenum], w_m;
	long_uint operanda[corenum], operandb[corenum], operandaout[corenum], operandbout[corenum];
	long_uint operanda_[corenum], operandb_[corenum], operandaout_[corenum], operandbout_[corenum];


ntt_stage1:
	for (uint32_t i = 1; i <= L_BRAMSIZE; i++)
	{
		w[0] = ONE; // no need to use w[1], in ntt_stage1
		w_m = OMEGA[i]; // w_m is 2^i-th root of unity
		stepsize = 1 << (i-1);
		stepnum = bramsize;

		MEa = 0;
		MEb = stepsize;

	ntt_stage1_inner:
		for (uint32_t j = 0; j < stepnum; j+=2)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = in inter false

			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				operanda[m] = in[m][MEa];
				operandb[m] = in[m][MEb];
				operanda_[m] = in[m + corenum][MEa];
				operandb_[m] = in[m + corenum][MEb];
			}

			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				// ntt core
				ntt_core_large_bw(w[0],p, operanda[m], operandb[m], operandaout[m], operandbout[m]);
				ntt_core_large_bw(w[0],p, operanda_[m], operandb_[m], operandaout_[m], operandbout_[m]);
			}

			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				in[m][MEa] = operandaout[m];
				in[m][MEb] = operandbout[m];
				in[m + corenum][MEa] = operandaout_[m];
				in[m + corenum][MEb] = operandbout_[m];
			}

			if(((MEa + 1) & stepsize)){
				MEa += stepsize + 1;
				MEb += stepsize + 1;
				w[0] = ONE;
			}
			else{
				MEa += 1;
				MEb += 1;
				w[0] = MUL_MOD(w[0], w_m,p);  // compute w for next ntt_core
			}
		}
	}
ntt_stage2:for(uint32_t i = L_BRAMSIZE + 1; i <= logn; i++) {
		//precompute w
		switch(i){
			case 8:
				w[0] = ONE;
				w[1] = ONE;
				w[2] = ONE;
				w[3] = ONE;
				w[4] = ONE;
				w[5] = ONE;
				w[6] = ONE;
				w[7] = ONE;
				w[8] = ONE;
				w[9] = ONE;
				w[10] = ONE;
				w[11] = ONE;
				w[12] = ONE;
				w[13] = ONE;
				w[14] = ONE;
				w[15] = ONE;
				break;
			case 9:
				w[0] = ONE;
				w[1] = OMEGA[2];
				w[1] =OMEGA[2];
				w[2] = ONE;
				w[3] = OMEGA[2];
				w[4] = ONE;
				w[5] = OMEGA[2];
				w[6] = ONE;
				w[7] = OMEGA[2];
				w[8] = ONE;
				w[9] = OMEGA[2];
				w[10] = ONE;
				w[11] = OMEGA[2];
				w[12] = ONE;
				w[13] = OMEGA[2];
				w[14] = ONE;
				w[15] = OMEGA[2];
				break;
			case 10:
				w[0] = ONE;
				w[1] = OMEGA[2];
				w[2] = ONE;
				w[3] = OMEGA[2];
				w[4] = ONE;
				w[5] = OMEGA[2];
				w[6] = ONE;
				w[7] = OMEGA[2];
				w[8] = ONE;
				w[9] = OMEGA[2];
				w[10] = ONE;
				w[11] = OMEGA[2];
				w[12] = ONE;
				w[13] = OMEGA[2];
				w[14] = ONE;
				w[15] = OMEGA[2];
				break;
			case 11:
				w[0] = ONE;
				w[1] = OMEGA[2];
				w[2] = ONE;
				w[3] = OMEGA[2];
				w[4] = ONE;
				w[5] = OMEGA[2];
				w[6] = ONE;
				w[7] = OMEGA[2];
				w[8] = ONE;
				w[9] = OMEGA[2];
				w[10] = ONE;
				w[11] = OMEGA[2];
				w[12] = ONE;
				w[13] = OMEGA[2];
				w[14] = ONE;
				w[15] = OMEGA[2];
				break;
			case 12:
				w[0] = ONE;
				w[1] = OMEGA[2];
				w[2] = ONE;
				w[3] = OMEGA[2];
				w[4] = ONE;
				w[5] = OMEGA[2];
				w[6] = ONE;
				w[7] = OMEGA[2];
				w[8] = ONE;
				w[9] = OMEGA[2];
				w[10] = ONE;
				w[11] = OMEGA[2];
				w[12] = ONE;
				w[13] = OMEGA[2];
				w[14] = ONE;
				w[15] = OMEGA[2];
		}
		w_m = OMEGA[i]; // w_m is 2^i-th root of unity
		stepsize = 1 << (i-1);
		stepnum = bramsize;
		MEa = 0;
		MEb = 0;
		ntt_stage2_inner:for(uint32_t j = 0; j < stepnum; j++){
#pragma HLS PIPELINE II = 1
#pragma HLS DEPENDENCE variable = in inter false
			switch(i){
				case 8:
					operanda[0] = in[0][MEa];
					operandb[0] = in[1][MEb];
					operanda[1] = in[2][MEa];
					operandb[1] = in[3][MEb];
					operanda[2] = in[4][MEa];
					operandb[2] = in[5][MEb];
					operanda[3] = in[6][MEa];
					operandb[3] = in[7][MEb];
					operanda[4] = in[8][MEa];
					operandb[4] = in[9][MEb];
					operanda[5] = in[10][MEa];
					operandb[5] = in[11][MEb];
					operanda[6] = in[12][MEa];
					operandb[6] = in[13][MEb];
					operanda[7] = in[14][MEa];
					operandb[7] = in[15][MEb];
					operanda[8] = in[16][MEa];
					operandb[8] = in[17][MEb];
					operanda[9] = in[18][MEa];
					operandb[9] = in[19][MEb];
					operanda[10] = in[20][MEa];
					operandb[10] = in[21][MEb];
					operanda[11] = in[22][MEa];
					operandb[11] = in[23][MEb];
					operanda[12] = in[24][MEa];
					operandb[12] = in[25][MEb];
					operanda[13] = in[26][MEa];
					operandb[13] = in[27][MEb];
					operanda[14] = in[28][MEa];
					operandb[14] = in[29][MEb];
					operanda[15] = in[30][MEa];
					operandb[15] = in[31][MEb];

					break;

				case 9:
					operanda[0] = in[0][MEa];
					operandb[0] = in[2][MEb];
					operanda[1] = in[1][MEa];
					operandb[1] = in[3][MEb];
					operanda[2] = in[4][MEa];
					operandb[2] = in[6][MEb];
					operanda[3] = in[5][MEa];
					operandb[3] = in[7][MEb];
					operanda[4] = in[8][MEa];
					operandb[4] = in[10][MEb];
					operanda[5] = in[9][MEa];
					operandb[5] = in[11][MEb];
					operanda[6] = in[12][MEa];
					operandb[6] = in[14][MEb];
					operanda[7] = in[13][MEa];
					operandb[7] = in[15][MEb];
					operanda[8] = in[16][MEa];
					operandb[8] = in[18][MEb];
					operanda[9] = in[17][MEa];
					operandb[9] = in[19][MEb];
					operanda[10] = in[20][MEa];
					operandb[10] = in[22][MEb];
					operanda[11] = in[21][MEa];
					operandb[11] = in[23][MEb];
					operanda[12] = in[24][MEa];
					operandb[12] = in[26][MEb];
					operanda[13] = in[25][MEa];
					operandb[13] = in[27][MEb];
					operanda[14] = in[28][MEa];
					operandb[14] = in[30][MEb];
					operanda[15] = in[29][MEa];
					operandb[15] = in[31][MEb];
					break;
				case 10:
					operanda[0] = in[0][MEa];
					operandb[0] = in[4][MEb];
					operanda[1] = in[1][MEa];
					operandb[1] = in[5][MEb];
					operanda[2] = in[2][MEa];
					operandb[2] = in[6][MEb];
					operanda[3] = in[3][MEa];
					operandb[3] = in[7][MEb];
					operanda[4] = in[8][MEa];
					operandb[4] = in[12][MEb];
					operanda[5] = in[9][MEa];
					operandb[5] = in[13][MEb];
					operanda[6] = in[10][MEa];
					operandb[6] = in[14][MEb];
					operanda[7] = in[11][MEa];
					operandb[7] = in[15][MEb];
					operanda[8] = in[16][MEa];
					operandb[8] = in[20][MEb];
					operanda[9] = in[17][MEa];
					operandb[9] = in[21][MEb];
					operanda[10] = in[18][MEa];
					operandb[10] = in[22][MEb];
					operanda[11] = in[19][MEa];
					operandb[11] = in[23][MEb];
					operanda[12] = in[24][MEa];
					operandb[12] = in[28][MEb];
					operanda[13] = in[25][MEa];
					operandb[13] = in[29][MEb];
					operanda[14] = in[26][MEa];
					operandb[14] = in[30][MEb];
					operanda[15] = in[27][MEa];
					operandb[15] = in[31][MEb];

					break;
				case 11:
					operanda[0] = in[0][MEa];
					operandb[0] = in[8][MEb];
					operanda[1] = in[1][MEa];
					operandb[1] = in[9][MEb];
					operanda[2] = in[2][MEa];
					operandb[2] = in[10][MEb];
					operanda[3] = in[3][MEa];
					operandb[3] = in[11][MEb];
					operanda[4] = in[4][MEa];
					operandb[4] = in[12][MEb];
					operanda[5] = in[5][MEa];
					operandb[5] = in[13][MEb];
					operanda[6] = in[6][MEa];
					operandb[6] = in[14][MEb];
					operanda[7] = in[7][MEa];
					operandb[7] = in[15][MEb];
					operanda[8] = in[16][MEa];
					operandb[8] = in[24][MEb];
					operanda[9] = in[17][MEa];
					operandb[9] = in[25][MEb];
					operanda[10] = in[18][MEa];
					operandb[10] = in[26][MEb];
					operanda[11] = in[19][MEa];
					operandb[11] = in[27][MEb];
					operanda[12] = in[20][MEa];
					operandb[12] = in[28][MEb];
					operanda[13] = in[21][MEa];
					operandb[13] = in[29][MEb];
					operanda[14] = in[22][MEa];
					operandb[14] = in[30][MEb];
					operanda[15] = in[23][MEa];
					operandb[15] = in[31][MEb];

					break;
				case 12:
					operanda[0] = in[0][MEa];
					operandb[0] = in[16][MEb];
					operanda[1] = in[1][MEa];
					operandb[1] = in[17][MEb];
					operanda[2] = in[2][MEa];
					operandb[2] = in[18][MEb];
					operanda[3] = in[3][MEa];
					operandb[3] = in[19][MEb];
					operanda[4] = in[4][MEa];
					operandb[4] = in[20][MEb];
					operanda[5] = in[5][MEa];
					operandb[5] = in[21][MEb];
					operanda[6] = in[6][MEa];
					operandb[6] = in[22][MEb];
					operanda[7] = in[7][MEa];
					operandb[7] = in[23][MEb];
					operanda[8] = in[8][MEa];
					operandb[8] = in[24][MEb];
					operanda[9] = in[9][MEa];
					operandb[9] = in[25][MEb];
					operanda[10] = in[10][MEa];
					operandb[10] = in[26][MEb];
					operanda[11] = in[11][MEa];
					operandb[11] = in[27][MEb];
					operanda[12] = in[12][MEa];
					operandb[12] = in[28][MEb];
					operanda[13] = in[13][MEa];
					operandb[13] = in[29][MEb];
					operanda[14] = in[14][MEa];
					operandb[14] = in[30][MEb];
					operanda[15] = in[15][MEa];
					operandb[15] = in[31][MEb];
			}

			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				// ntt core
				ntt_core_large_bw(w[m],p, operanda[m], operandb[m], operandaout[m], operandbout[m]);
			}


			switch(i){
				case 8:
					in[0][MEa] = operandaout[0];
					in[1][MEb] = operandbout[0];
					in[2][MEa] = operandaout[1];
					in[3][MEb] = operandbout[1];
					in[4][MEa] = operandaout[2];
					in[5][MEb] = operandbout[2];
					in[6][MEa] = operandaout[3];
					in[7][MEb] = operandbout[3];
					in[8][MEa] = operandaout[4];
					in[9][MEb] = operandbout[4];
					in[10][MEa] = operandaout[5];
					in[11][MEb] = operandbout[5];
					in[12][MEa] = operandaout[6];
					in[13][MEb] = operandbout[6];
					in[14][MEa] = operandaout[7];
					in[15][MEb] = operandbout[7];
					in[16][MEa] = operandaout[8];
					in[17][MEb] = operandbout[8];
					in[18][MEa] = operandaout[9];
					in[19][MEb] = operandbout[9];
					in[20][MEa] = operandaout[10];
					in[21][MEb] = operandbout[10];
					in[22][MEa] = operandaout[11];
					in[23][MEb] = operandbout[11];
					in[24][MEa] = operandaout[12];
					in[25][MEb] = operandbout[12];
					in[26][MEa] = operandaout[13];
					in[27][MEb] = operandbout[13];
					in[28][MEa] = operandaout[14];
					in[29][MEb] = operandbout[14];
					in[30][MEa] = operandaout[15];
					in[31][MEb] = operandbout[15];
				break;
				case 9:
					in[0][MEa] = operandaout[0];
					in[2][MEb] = operandbout[0];
					in[1][MEa] = operandaout[1];
					in[3][MEb] = operandbout[1];
					in[4][MEa] = operandaout[2];
					in[6][MEb] = operandbout[2];
					in[5][MEa] = operandaout[3];
					in[7][MEb] = operandbout[3];
					in[8][MEa] = operandaout[4];
					in[10][MEb] = operandbout[4];
					in[9][MEa] = operandaout[5];
					in[11][MEb] = operandbout[5];
					in[12][MEa] = operandaout[6];
					in[14][MEb] = operandbout[6];
					in[13][MEa] = operandaout[7];
					in[15][MEb] = operandbout[7];
					in[16][MEa] = operandaout[8];
					in[18][MEb] = operandbout[8];
					in[17][MEa] = operandaout[9];
					in[19][MEb] = operandbout[9];
					in[20][MEa] = operandaout[10];
					in[22][MEb] = operandbout[10];
					in[21][MEa] = operandaout[11];
					in[23][MEb] = operandbout[11];
					in[24][MEa] = operandaout[12];
					in[26][MEb] = operandbout[12];
					in[25][MEa] = operandaout[13];
					in[27][MEb] = operandbout[13];
					in[28][MEa] = operandaout[14];
					in[30][MEb] = operandbout[14];
					in[29][MEa] = operandaout[15];
					in[31][MEb] = operandbout[15];
					break;
				case 10:
					in[0][MEa] = operandaout[0];
					in[4][MEb] = operandbout[0];
					in[1][MEa] = operandaout[1];
					in[5][MEb] = operandbout[1];
					in[2][MEa] = operandaout[2];
					in[6][MEb] = operandbout[2];
					in[3][MEa] = operandaout[3];
					in[7][MEb] = operandbout[3];
					in[8][MEa] = operandaout[4];
					in[12][MEb] = operandbout[4];
					in[9][MEa] = operandaout[5];
					in[13][MEb] = operandbout[5];
					in[10][MEa] = operandaout[6];
					in[14][MEb] = operandbout[6];
					in[11][MEa] = operandaout[7];
					in[15][MEb] = operandbout[7];
					in[16][MEa] = operandaout[8];
					in[20][MEb] = operandbout[8];
					in[17][MEa] = operandaout[9];
					in[21][MEb] = operandbout[9];
					in[18][MEa] = operandaout[10];
					in[22][MEb] = operandbout[10];
					in[19][MEa] = operandaout[11];
					in[23][MEb] = operandbout[11];
					in[24][MEa] = operandaout[12];
					in[28][MEb] = operandbout[12];
					in[25][MEa] = operandaout[13];
					in[29][MEb] = operandbout[13];
					in[26][MEa] = operandaout[14];
					in[30][MEb] = operandbout[14];
					in[27][MEa] = operandaout[15];
					in[31][MEb] = operandbout[15];

					break;
				case 11:
					in[0][MEa] = operandaout[0];
					in[8][MEb] = operandbout[0];
					in[1][MEa] = operandaout[1];
					in[9][MEb] = operandbout[1];
					in[2][MEa] = operandaout[2];
					in[10][MEb] = operandbout[2];
					in[3][MEa] = operandaout[3];
					in[11][MEb] = operandbout[3];
					in[4][MEa] = operandaout[4];
					in[12][MEb] = operandbout[4];
					in[5][MEa] = operandaout[5];
					in[13][MEb] = operandbout[5];
					in[6][MEa] = operandaout[6];
					in[14][MEb] = operandbout[6];
					in[7][MEa] = operandaout[7];
					in[15][MEb] = operandbout[7];
					in[16][MEa] = operandaout[8];
					in[24][MEb] = operandbout[8];
					in[17][MEa] = operandaout[9];
					in[25][MEb] = operandbout[9];
					in[18][MEa] = operandaout[10];
					in[26][MEb] = operandbout[10];
					in[19][MEa] = operandaout[11];
					in[27][MEb] = operandbout[11];
					in[20][MEa] = operandaout[12];
					in[28][MEb] = operandbout[12];
					in[21][MEa] = operandaout[13];
					in[29][MEb] = operandbout[13];
					in[22][MEa] = operandaout[14];
					in[30][MEb] = operandbout[14];
					in[23][MEa] = operandaout[15];
					in[31][MEb] = operandbout[15];

					break;
				case 12:
					in[0][MEa] = operandaout[0];
					in[16][MEb] = operandbout[0];
					in[1][MEa] = operandaout[1];
					in[17][MEb] = operandbout[1];
					in[2][MEa] = operandaout[2];
					in[18][MEb] = operandbout[2];
					in[3][MEa] = operandaout[3];
					in[19][MEb] = operandbout[3];
					in[4][MEa] = operandaout[4];
					in[20][MEb] = operandbout[4];
					in[5][MEa] = operandaout[5];
					in[21][MEb] = operandbout[5];
					in[6][MEa] = operandaout[6];
					in[22][MEb] = operandbout[6];
					in[7][MEa] = operandaout[7];
					in[23][MEb] = operandbout[7];
					in[8][MEa] = operandaout[8];
					in[24][MEb] = operandbout[8];
					in[9][MEa] = operandaout[9];
					in[25][MEb] = operandbout[9];
					in[10][MEa] = operandaout[10];
					in[26][MEb] = operandbout[10];
					in[11][MEa] = operandaout[11];
					in[27][MEb] = operandbout[11];
					in[12][MEa] = operandaout[12];
					in[28][MEb] = operandbout[12];
					in[13][MEa] = operandaout[13];
					in[29][MEb] = operandbout[13];
					in[14][MEa] = operandaout[14];
					in[30][MEb] = operandbout[14];
					in[15][MEa] = operandaout[15];
					in[31][MEb] = operandbout[15];


			MEa += 1;
			MEb += 1;

			for(int m = 0; m < corenum; m++){
#pragma HLS UNROLL
				w[m] = MUL_MOD(w[m], w_m,p);  // compute w for next ntt_core
			}


		}
	}
}
}
