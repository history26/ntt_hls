#include "ntt_v2.hpp"


void ntt_4core_new(long_uint operand[bramnum][bramsize],
				   long_uint rp[RPBRAMNUM][RPBRAMSIZE],
				   long_uint modulus)
{
#pragma HLS INLINE off

	long_uint operanda[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	long_uint operandb[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	long_uint operanda_[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	long_uint operandb_[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	long_uint operanda_out[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_out complete dim = 0
	long_uint operandb_out[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_out complete dim = 0
	long_uint operanda_2out[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_2out complete dim = 0
	long_uint operandb_2out[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_2out complete dim = 0

	long_uint two_times_modulus = modulus << 1;
	long_uint one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

ntt_stage1:
	for (uint32_t i = 0; i < stage1_max; i++)
	{
		uint32_t stepsize = step_num >> (i + 1);

		uint32_t temp1 = 1 << i;
		uint32_t temp2 = stage1_max - i;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;

		MEa = 0;
		MEb = stepsize;

	ntt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			IDXTYPE rp_idx0 = temp1 + (j >> temp2);
			IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
			IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

			//			rp_idx0_i = rp_idx0 / RPBRAMNUM;
			//			rp_idx0_j = rp_idx0 % RPBRAMSIZE;
			//			rp_idx00_i = rp_idx00 / RPBRAMNUM;
			//			rp_idx00_j = rp_idx00 % RPBRAMSIZE;

			rp_idx0_i = rp_idx0(L_RPBRAMNUM - 1, 0);
			rp_idx0_j = rp_idx0 >> L_RPBRAMNUM;
			rp_idx00_i = rp_idx00(L_RPBRAMNUM - 1, 0);
			rp_idx00_j = rp_idx00 >> L_RPBRAMNUM;

			long_uint w0 = rp[rp_idx0_i][rp_idx0_j];
			// long_uint sw0 = srp[rp_idx0_i][rp_idx0_j];
			long_uint w00 = rp[rp_idx00_i][rp_idx00_j];
			// long_uint sw00 = srp[rp_idx00_i][rp_idx00_j];

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				operanda[m] = operand[m][MEa];
				operandb[m] = operand[m][MEb];
				operanda_[m] = operand[m + corenum][MEa];
				operandb_[m] = operand[m + corenum][MEb];
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL

				// ntt_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m],
				// 		 operanda_out[m], operandb_out[m]);
                ntt_core_large_bw(w0, modulus, operanda[m],operandb[m],
						operanda_out[m],operandb_out[m]);

				// ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[m], operandb_[m],
				// 		 operanda_2out[m], operandb_2out[m]);
                ntt_core_large_bw(w00, modulus, operanda_[m],operandb_[m],
						operanda_2out[m],operandb_2out[m]);
			}
			j++;

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				operand[m][MEa] = operanda_out[m];
				operand[m][MEb] = operandb_out[m];
				operand[m + corenum][MEa] = operanda_2out[m];
				operand[m + corenum][MEb] = operandb_2out[m];
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}

ntt_stage2_0:
	for (uint32_t i = stage1_max; i < stage1_max + 1; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> temp4;
		uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << i) + (m >> temp4);

			// cout << "arrayWindex["<< m << "] = "<<arrayWindex[m] << endl;
		}

	ntt_stage2_inner_0:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			long_uint w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
// 			long_uint sw[corenum];
// #pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;
			// cout << "temp_adder = " << temp_adder << endl;

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL

				IDXTYPE x = (arrayWindex[0] + temp_adder);
				IDXTYPE rp_idx0_i, rp_idx0_j;
				rp_idx0_i = x(L_RPBRAMNUM - 1, 0);
				rp_idx0_j = x >> L_RPBRAMNUM;
				// cout << "x = " << x << endl;

				w[m] = rp[rp_idx0_i][rp_idx0_j];
				// sw[m] = srp[rp_idx0_i][rp_idx0_j];
			}

			operanda[0] = operand[0][MEa];
			operandb[0] = operand[4][MEb];
			operanda[1] = operand[1][MEa];
			operandb[1] = operand[5][MEb];
			operanda[2] = operand[2][MEa];
			operandb[2] = operand[6][MEb];
			operanda[3] = operand[3][MEa];
			operandb[3] = operand[7][MEb];

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL

				// cout << "w[m]" << w[m] << endl;
				// cout << "sw[m]" << sw[m] << endl;

				// cout << "in: operanda[m]" << operanda[m] << endl;
				// cout << "in: operandb[m]" << operandb[m] << endl;
				// ntt_core(w[m], sw[m], modulus, two_times_modulus, operanda[m], operandb[m], operanda_[m], operandb_[m]);
                
                ntt_core_large_bw(w[m], modulus, operanda[m],operandb[m],
						operanda_[m],operandb_[m]);

				// cout << "out: operanda[m]" << operanda_[m] << endl;
				// cout << "out: operandb[m]" << operandb_[m] << endl;
			}

			operand[0][MEa] = operanda_[0];
			operand[4][MEb] = operandb_[0];
			operand[1][MEa] = operanda_[1];
			operand[5][MEb] = operandb_[1];
			operand[2][MEa] = operanda_[2];
			operand[6][MEb] = operandb_[2];
			operand[3][MEa] = operanda_[3];
			operand[7][MEb] = operandb_[3];

			MEa += 1;
			MEb += 1;
		}
	}

ntt_stage2_1:
	for (uint32_t i = stage1_max + 1; i < stage1_max + 2; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> temp4;
		uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << i) + (m >> temp4);
		}

	ntt_stage2_inner_1:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			long_uint w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
// 			long_uint sw[corenum];
// #pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			IDXTYPE x = (arrayWindex[0] + temp_adder);
			IDXTYPE rp_idx0_i, rp_idx0_j;
			rp_idx0_i = x(L_RPBRAMNUM - 1, 0);
			rp_idx0_j = x >> L_RPBRAMNUM;

			long_uint rp_t1, srp_t1;
			long_uint rp_t2, srp_t2;

			rp_t1 = rp[rp_idx0_i][rp_idx0_j];
			// srp_t1 = srp[rp_idx0_i][rp_idx0_j];

			rp_t2 = rp[rp_idx0_i + 1][rp_idx0_j];
			// srp_t2 = srp[rp_idx0_i + 1][rp_idx0_j];

			w[0] = rp_t1;
			// sw[0] = srp_t1;
			w[1] = rp_t1;
			// sw[1] = srp_t1;
			w[2] = rp_t2;
			// sw[2] = srp_t2;
			w[3] = rp_t2;
			// sw[3] = srp_t2;

			operanda[0] = operand[0][MEa];
			operandb[0] = operand[2][MEb];
			operanda[1] = operand[1][MEa];
			operandb[1] = operand[3][MEb];
			operanda[2] = operand[4][MEa];
			operandb[2] = operand[6][MEb];
			operanda[3] = operand[5][MEa];
			operandb[3] = operand[7][MEb];

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL

				// ntt_core(w[m], sw[m], modulus, two_times_modulus, operanda[m], operandb[m], operanda_[m], operandb_[m]);
				ntt_core_large_bw(w[m], modulus, operanda[m],operandb[m],
						operanda_[m],operandb_[m]);
			}

			operand[0][MEa] = operanda_[0];
			operand[2][MEb] = operandb_[0];
			operand[1][MEa] = operanda_[1];
			operand[3][MEb] = operandb_[1];
			operand[4][MEa] = operanda_[2];
			operand[6][MEb] = operandb_[2];
			operand[5][MEa] = operanda_[3];
			operand[7][MEb] = operandb_[3];

			MEa += 1;
			MEb += 1;
		}
	}

ntt_stage2_2:
	for (uint32_t i = stage1_max + 2; i < stage_num; i++)
	{

		uint32_t temp4 = stage_num - 1 - i;
		uint32_t temp_par = corenum >> temp4;
		uint32_t temp5 = i - stage1_max;

		MEa = 0;
		MEb = 0;

		for (int m = 0; m < corenum; m++)
		{
			arrayWindex[m] = (1 << i) + (m >> temp4);
		}

	ntt_stage2_inner_2:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable = operand inter false
#pragma HLS DEPENDENCE variable = operanda inter false
#pragma HLS DEPENDENCE variable = operandb inter false

			long_uint w[corenum];
#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
// 			long_uint sw[corenum];
// #pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

			uint32_t temp_adder = temp_par * j;

			IDXTYPE x = (arrayWindex[0] + temp_adder);
			IDXTYPE rp_idx0_j;
			rp_idx0_j = x >> L_RPBRAMNUM;

			w[0] = rp[0][rp_idx0_j];
			// sw[0] = srp[0][rp_idx0_j];
			w[1] = rp[1][rp_idx0_j];
			// sw[1] = srp[1][rp_idx0_j];
			w[2] = rp[2][rp_idx0_j];
			// sw[2] = srp[2][rp_idx0_j];
			w[3] = rp[3][rp_idx0_j];
			// sw[3] = srp[3][rp_idx0_j];

			operanda[0] = operand[0][MEa];
			operandb[0] = operand[1][MEb];
			operanda[1] = operand[2][MEa];
			operandb[1] = operand[3][MEb];
			operanda[2] = operand[4][MEa];
			operandb[2] = operand[5][MEb];
			operanda[3] = operand[6][MEa];
			operandb[3] = operand[7][MEb];

			for (uint32_t m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				// ntt_core(w[m], sw[m], modulus, two_times_modulus,
				// 		 operanda[m], operandb[m], operanda_[m], operandb_[m]);
				ntt_core_large_bw(w[m], modulus, operanda[m],operandb[m],
						operanda_[m],operandb_[m]);
			}

			operand[0][MEa] = operanda_[0];
			operand[1][MEb] = operandb_[0];
			operand[2][MEa] = operanda_[1];
			operand[3][MEb] = operandb_[1];
			operand[4][MEa] = operanda_[2];
			operand[5][MEb] = operandb_[2];
			operand[6][MEa] = operanda_[3];
			operand[7][MEb] = operandb_[3];

			MEa += 1;
			MEb += 1;
		}
	}
}


// change the stage_2 loop


void ntt_4core_new_v2(long_uint operand[bramnum][bramsize],
				   long_uint rp[RPBRAMNUM][RPBRAMSIZE],
				   long_uint modulus)
{
#pragma HLS INLINE off

	long_uint operanda[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda complete dim = 0
	long_uint operandb[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb complete dim = 0
	long_uint operanda_[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_ complete dim = 0
	long_uint operandb_[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_ complete dim = 0

	long_uint operanda_out[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_out complete dim = 0
	long_uint operandb_out[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_out complete dim = 0
	long_uint operanda_2out[corenum];
#pragma HLS ARRAY_PARTITION variable = operanda_2out complete dim = 0
	long_uint operandb_2out[corenum];
#pragma HLS ARRAY_PARTITION variable = operandb_2out complete dim = 0

	long_uint two_times_modulus = modulus << 1;
	long_uint one_u64 = 1;
	uint32_t step_num = bramsize;
	uint32_t stage_num = STAGENUM;
	uint32_t stage1_max = stagemax;

	uint32_t MEa;
	uint32_t MEb;

	uint32_t arrayWindex[corenum];
#pragma HLS ARRAY_PARTITION variable = arrayWindex complete dim = 0

ntt_stage1:
	for (uint32_t i = 0; i < stage1_max; i++)
	{
		uint32_t stepsize = step_num >> (i + 1);

		uint32_t temp1 = 1 << i;
		uint32_t temp2 = stage1_max - i;
		uint32_t temp3 = (one_u64 << temp2) - one_u64;

		MEa = 0;
		MEb = stepsize;

	ntt_stage1_inner:
		for (uint32_t j = 0; j < step_num; j++)
		{
#pragma HLS PIPELINE II = 2
#pragma HLS DEPENDENCE variable = operand inter false

			IDXTYPE rp_idx0 = temp1 + (j >> temp2);
			IDXTYPE rp_idx00 = temp1 + ((j + 1) >> temp2);
			IDXTYPE rp_idx0_i, rp_idx0_j, rp_idx00_i, rp_idx00_j;

			//			rp_idx0_i = rp_idx0 / RPBRAMNUM;
			//			rp_idx0_j = rp_idx0 % RPBRAMSIZE;
			//			rp_idx00_i = rp_idx00 / RPBRAMNUM;
			//			rp_idx00_j = rp_idx00 % RPBRAMSIZE;

			rp_idx0_i = rp_idx0(L_RPBRAMNUM - 1, 0);
			rp_idx0_j = rp_idx0 >> L_RPBRAMNUM;
			rp_idx00_i = rp_idx00(L_RPBRAMNUM - 1, 0);
			rp_idx00_j = rp_idx00 >> L_RPBRAMNUM;

			long_uint w0 = rp[rp_idx0_i][rp_idx0_j];
			// long_uint sw0 = srp[rp_idx0_i][rp_idx0_j];
			long_uint w00 = rp[rp_idx00_i][rp_idx00_j];
			// long_uint sw00 = srp[rp_idx00_i][rp_idx00_j];

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				operanda[m] = operand[m][MEa];
				operandb[m] = operand[m][MEb];
				operanda_[m] = operand[m + corenum][MEa];
				operandb_[m] = operand[m + corenum][MEb];
			}

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL

				// ntt_core(w0, sw0, modulus, two_times_modulus, operanda[m], operandb[m],
				// 		 operanda_out[m], operandb_out[m]);
				ntt_core_large_bw(w0, modulus, operanda[m],operandb[m],
						operanda_out[m],operandb_out[m]);

				// ntt_core(w00, sw00, modulus, two_times_modulus, operanda_[m], operandb_[m],
				// 		 operanda_2out[m], operandb_2out[m]);
				ntt_core_large_bw(w00, modulus, operanda_[m],operandb_[m],
						operanda_2out[m],operandb_2out[m]);
			}
			j++;

			for (int m = 0; m < corenum; m++)
			{
#pragma HLS UNROLL
				operand[m][MEa] = operanda_out[m];
				operand[m][MEb] = operandb_out[m];
				operand[m + corenum][MEa] = operanda_2out[m];
				operand[m + corenum][MEb] = operandb_2out[m];
			}

			if (((j + 1) & temp3) == 0)
			{
				MEa += stepsize + 1;
				MEb += stepsize + 1;
			}
			else
			{
				MEa += 1;
				MEb += 1;
			}
		}
	}

ntt_stage2:
	for(uint32_t i=stage1_max;i<stage_num;i++)
	{

	#pragma HLS LOOP_FLATTEN off
			uint32_t temp4 = stage_num - 1 - i;
			uint32_t temp_par = corenum >> temp4;
			uint32_t temp5 = i - stage1_max;

			MEa = 0;
			MEb = 0;

			for (int m = 0; m < corenum; m++)
			{
				arrayWindex[m] = (1 << i) + (m >> temp4);
			}

		ntt_stage2_inner:
			for (uint32_t j = 0; j < step_num; j++)
			{
			#pragma HLS PIPELINE
			#pragma HLS DEPENDENCE variable = operand inter false
			#pragma HLS DEPENDENCE variable = operanda inter false
			#pragma HLS DEPENDENCE variable = operandb inter false

				long_uint w[corenum];
			#pragma HLS ARRAY_PARTITION variable = w complete dim = 1
			// 			long_uint sw[corenum];
			// #pragma HLS ARRAY_PARTITION variable = sw complete dim = 1

				uint32_t temp_adder = temp_par * j;

				for (uint32_t m = 0; m < corenum; m++)
				{
	#pragma HLS UNROLL

					IDXTYPE x = (arrayWindex[0] + temp_adder);
					IDXTYPE rp_idx0_i, rp_idx0_j;
					rp_idx0_i = x(L_RPBRAMNUM - 1, 0);
					rp_idx0_j = x >> L_RPBRAMNUM;
					// cout << "x = " << x << endl;

					w[m] = rp[rp_idx0_i][rp_idx0_j];
					// sw[m] = srp[rp_idx0_i][rp_idx0_j];
				}

				switch (i)
				{
				case (stagemax):
					operanda[0] = operand[0][MEa];
					operandb[0] = operand[4][MEb];
					operanda[1] = operand[1][MEa];
					operandb[1] = operand[5][MEb];
					operanda[2] = operand[2][MEa];
					operandb[2] = operand[6][MEb];
					operanda[3] = operand[3][MEa];
					operandb[3] = operand[7][MEb];
					break;

				case (stagemax+1):
					operanda[0] = operand[0][MEa];
					operandb[0] = operand[2][MEb];
					operanda[1] = operand[1][MEa];
					operandb[1] = operand[3][MEb];
					operanda[2] = operand[4][MEa];
					operandb[2] = operand[6][MEb];
					operanda[3] = operand[5][MEa];
					operandb[3] = operand[7][MEb];
					break;
				
				case (stagemax+2):
					operanda[0] = operand[0][MEa];
					operandb[0] = operand[1][MEb];
					operanda[1] = operand[2][MEa];
					operandb[1] = operand[3][MEb];
					operanda[2] = operand[4][MEa];
					operandb[2] = operand[5][MEb];
					operanda[3] = operand[6][MEa];
					operandb[3] = operand[7][MEb];
				default:
					break;
				}

				for (uint32_t m = 0; m < corenum; m++)
				{
	#pragma HLS UNROLL

					// cout << "w[m]" << w[m] << endl;
					// cout << "sw[m]" << sw[m] << endl;

					// cout << "in: operanda[m]" << operanda[m] << endl;
					// cout << "in: operandb[m]" << operandb[m] << endl;
					// ntt_core(w[m], sw[m], modulus, two_times_modulus, operanda[m], operandb[m], operanda_[m], operandb_[m]);
					
					ntt_core_large_bw(w[m], modulus, operanda[m],operandb[m],
							operanda_[m],operandb_[m]);

					// cout << "out: operanda[m]" << operanda_[m] << endl;
					// cout << "out: operandb[m]" << operandb_[m] << endl;
				}

				
				switch (i)
				{
				case (stagemax):
					operand[0][MEa] = operanda_[0];
					operand[4][MEb] = operandb_[0];
					operand[1][MEa] = operanda_[1];
					operand[5][MEb] = operandb_[1];
					operand[2][MEa] = operanda_[2];
					operand[6][MEb] = operandb_[2];
					operand[3][MEa] = operanda_[3];
					operand[7][MEb] = operandb_[3];
					break;

				case (stagemax+1):
					operand[0][MEa] = operanda_[0];
					operand[2][MEb] = operandb_[0];
					operand[1][MEa] = operanda_[1];
					operand[3][MEb] = operandb_[1];
					operand[4][MEa] = operanda_[2];
					operand[6][MEb] = operandb_[2];
					operand[5][MEa] = operanda_[3];
					operand[7][MEb] = operandb_[3];
					break;

				case(stagemax+2):
					operand[0][MEa] = operanda_[0];
					operand[1][MEb] = operandb_[0];
					operand[2][MEa] = operanda_[1];
					operand[3][MEb] = operandb_[1];
					operand[4][MEa] = operanda_[2];
					operand[5][MEb] = operandb_[2];
					operand[6][MEa] = operanda_[3];
					operand[7][MEb] = operandb_[3];
				
				default:
					break;
				}
				MEa += 1;
				MEb += 1;
			}
	}

}
