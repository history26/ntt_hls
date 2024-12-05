#ifndef _FOUR_STEP_NTT_FUNC_H_
#define _FOUR_STEP_NTT_FUNC_H_

#include "ntt.hpp"


void Four_step_NTT_2_24_Ncore_large_bw(in_type* in,in_type* out,in_type* temp1,in_type* temp2);
void Test(in_type* in,in_type* out);
//void Test_Uram(in_type* in,in_type* out);
void Test_Uram(in_type* in,in_type* out);

//I/O functions

//read the position-th row,then put data into res with bit-reverse
void read_row_reverse(in_type* in,long_uint res[bramnum][bramsize],uint32_t position);
//read the position-th column,then put data into res with bit-reverse
void read_column_reverse(in_type* in,long_uint res[bramnum][bramsize],uint32_t position);

void read_column_reverse_para2(in_type* in,long_uint res[2][bramnum][bramsize],uint32_t position);

void step_4(in_type* out,long_uint res[bramnum][bramsize],uint32_t position);

void test_para(in_type* in1,in_type* out1,in_type* in2,in_type* out2);

// template <unsigned corenum, unsigned bramnum, unsigned bramsize, unsigned stagemax>
// void ntt_1core(UDTYPE operand[bramnum][bramsize],

#endif
