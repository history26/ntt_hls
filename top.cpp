#include "top.h"


void dut(in_type* in,in_type* out,in_type* temp1,in_type* temp2){
#pragma HLS INTERFACE s_axilite register port = return
#pragma HLS INTERFACE m_axi depth = 16777216 port = in offset = slave bundle = gmem0
#pragma HLS INTERFACE m_axi depth = 16777216 port = out offset = slave bundle = gmem1
#pragma HLS INTERFACE m_axi depth = 16777216 port = temp1 offset = slave bundle = gmem2
#pragma HLS INTERFACE m_axi depth = 16777216 port = temp2 offset = slave bundle = gmem3
// #pragma HLS INTERFACE m_axi depth = 16777216 port = temp3 offset = slave bundle = gmem4
// #pragma HLS INTERFACE m_axi depth = 16777216 port = temp4 offset = slave bundle = gmem5
// #pragma HLS INTERFACE m_axi depth = 16777216 port = temp5 offset = slave bundle = gmem6
// #pragma HLS INTERFACE m_axi depth = 16777216 port = temp6 offset = slave bundle = gmem7
// #pragma HLS INTERFACE m_axi depth = 16777216 port = temp7 offset = slave bundle = gmem8
// #pragma HLS INTERFACE m_axi depth = 16777216 port = temp8 offset = slave bundle = gmem9
// #pragma HLS INTERFACE m_axi depth = 16777216 port = temp9 offset = slave bundle = gmem10
// #pragma HLS INTERFACE m_axi depth = 16777216 port = temp10 offset = slave bundle = gmem11



	// ntt_4core_new_v2(step1_in, rp,p);
	test_para(in,out,temp1,temp2);


	return;
}
