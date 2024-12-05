#include "ntt.hpp"
#include "ntt_v2.hpp"
#include "four_step_ntt.hpp"
// #include "ap_int.h"


#define in_bw 256
#define out_bw 512

//#define N 1024
//#define BASE_WIDTH 256
//#define in_WIDTH 512
//typedef ap_uint<in_WIDTH> in_type;
//typedef ap_uint<257> long_uint;
//typedef ap_uint<512> long_long_uint;
//typedef ap_int<512> long_int;


typedef ap_uint<out_bw> test_out;
typedef ap_uint<in_bw> test_in;

const ap_uint<in_bw> test_inv("0x04320FEDCBA987654320FEDCBA987654320FEDCBA987654320FEDCBA987654321",16);
const ap_uint<in_bw> test_omega[8]={
		ap_uint<in_bw>("0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb",16),
		ap_uint<in_bw>("0x04320FEDCBA987654320FEDCBA987654320FEDCBA987654320FEDCBA987654321",16),
		ap_uint<in_bw>("0x04320FEDCBA987654320FEDCBA987654320FEDCBA987654320FEDCBA987654321",16),
		ap_uint<in_bw>("0x04320FEDCBA987654320FEDCBA987654320FEDCBA987654320FEDCBA987654321",16),
		ap_uint<in_bw>("0x04320FEDCBA987654320FEDCBA987654320FEDCBA987654320FEDCBA987654321",16),
		ap_uint<in_bw>("0x04320FEDCBA987654320FEDCBA987654320FEDCBA987654320FEDCBA987654321",16),
		ap_uint<in_bw>("0x04320FEDCBA987654320FEDCBA987654320FEDCBA987654320FEDCBA987654324",16),
		ap_uint<in_bw>("0x04320FEDCBA987654320FEDCBA987654320FEDCBA987654320FEDCBA987654327",16)

};

void dut(in_type* in,in_type* out,in_type* temp1,in_type* temp2);
