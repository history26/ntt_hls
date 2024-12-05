#include "four_step_ntt.hpp"
using namespace std;

void read_row_reverse(in_type* in,long_uint res[bramnum][bramsize],uint32_t position){
	typedef ap_uint<12> INDEX;
	uint32_t start=position*C;
		for(int j=0;j<C;j++){
			INDEX rj=j;
			rj.reverse();
			INDEX idx_1,idx_2;
			idx_1=rj(11,12 - L_BRAMNUM);
			idx_2=rj(12 - L_BRAMNUM - 1,0);
			res[idx_1][idx_2]=in[start+j];
		}
}

void read_row_reverse_para2(in_type* in,long_uint res[2][bramnum][bramsize],uint32_t position){
	typedef ap_uint<12> INDEX;
	uint32_t trip_count=C/2;
	uint32_t start_1=position*trip_count;
	uint32_t start_2=start_1+trip_count;
	for(int j=0;j<trip_count;j++){
		#pragma HLS PIPELINE II = 2
		INDEX rj=2*j;
		rj.reverse();
		INDEX idx_1,idx_2;
		idx_1=rj(11,12 - L_BRAMNUM);
		idx_2=rj(12 - L_BRAMNUM - 1,0);

		//index of second element
		INDEX rj_2=2*j+1;
		rj_2.reverse();
		INDEX idx_1_2,idx_2_2;
		idx_1_2=rj_2(11,12 - L_BRAMNUM);
		idx_2_2=rj_2(12 - L_BRAMNUM - 1,0);
		//read element from first&second row ,2*2=4 elements
		in_type temp_1=in[start_1+j];
		in_type temp_2=in[start_2+j];

		res[0][idx_1][idx_2]=temp_1((2*BASE_WIDTH-1),BASE_WIDTH);
		res[0][idx_1_2][idx_2_2]=temp_1(BASE_WIDTH-1,0);
		res[1][idx_1][idx_2]=temp_2((2*BASE_WIDTH-1),BASE_WIDTH);
		res[1][idx_1_2][idx_2_2]=temp_2(BASE_WIDTH-1,0);
	}
}

void read_column_reverse(in_type* in,long_uint res[bramnum][bramsize],uint32_t position){
	typedef ap_uint<12> INDEX;
		for(int j=0;j<R;j++){
            #pragma HLS DEPENDENCE variable = res inter false
            #pragma HLS DEPENDENCE variable = res intra false
			INDEX rj=j;
			rj.reverse();
			INDEX idx_1,idx_2;
			idx_1=rj(11,12 - L_BRAMNUM);
			idx_2=rj(12 - L_BRAMNUM - 1,0);
			res[idx_1][idx_2]=in[j*C+position];
		}
}

void read_column_reverse_para2(in_type* in,long_uint res[2][bramnum][bramsize],uint32_t position){
	//Assume the matrix is stored in rows
	//para 2 means reading two column elements at same time

	typedef ap_uint<12> INDEX;
		for(int j=0;j<R;j++){
            #pragma HLS DEPENDENCE variable = res inter false
            #pragma HLS DEPENDENCE variable = res intra false

			INDEX rj=j;
			rj.reverse();
			INDEX idx_1,idx_2;
			idx_1=rj(11,12 - L_BRAMNUM);
			idx_2=rj(12 - L_BRAMNUM - 1,0);
            in_type temp=in[(j*C+position)/2];
			res[0][idx_1][idx_2]=temp((2*BASE_WIDTH-1),BASE_WIDTH);
            res[1][idx_1][idx_2]=temp(BASE_WIDTH-1,0);
		}
}


void step_2(in_type* out,long_uint res[bramnum][bramsize],uint32_t pre_w,uint32_t position){
		typedef ap_uint<12> INDEX;
		typedef ap_uint<32> BIG_INDEX;
		//step-2
		//Multiply each element by r^(i*j)
		long_uint w=ONE;
		step_2:for(BIG_INDEX j=0;j<R;j++){
			//i-th column ,j-th row
			INDEX idx_1,idx_2;
			idx_1=j(11,12 - L_BRAMNUM);
			idx_2=j(12 - L_BRAMNUM - 1,0);
			long_uint temp=MUL_MOD(res[idx_1][idx_2],w,p);
			out[j*C+position]=temp;
			w=MUL_MOD(w,pre_w,p);
		}
}

void step_2_para(in_type* out,long_uint res[2][bramnum][bramsize],uint32_t pre_w,uint32_t position){
		typedef ap_uint<12> INDEX;
		typedef ap_uint<32> BIG_INDEX;
		//step-2
		//
		long_uint w=ONE;
		long_uint w_after=ONE;
		long_uint w_1=OMEGA[24];
		long_uint pre_w_after_position=MUL_MOD(pre_w, w_1,p);
		step_2:for(BIG_INDEX j=0;j<R;j++){
			//Multiply each element by r^(i*j)
			//i-th column ,j-th row
			INDEX idx_1,idx_2;
			idx_1=j(11,12 - L_BRAMNUM);
			idx_2=j(12 - L_BRAMNUM - 1,0);
			in_type temp=MUL_MOD(res[0][idx_1][idx_2],w,p);
			temp=temp<<BASE_WIDTH;
			temp+=MUL_MOD(res[1][idx_1][idx_2],w_after,p);
			out[(j*C+position)/2]=temp;
			w=MUL_MOD(w,pre_w,p);
			w_after=MUL_MOD(w_after,pre_w_after_position,p);
		}
}

void step_4(in_type* out,long_uint res[bramnum][bramsize],uint32_t position){
		//transpose the matrix
		//write to position-th column
		typedef ap_uint<12> INDEX;
		typedef ap_uint<32> BIG_INDEX;
		step_4:for(BIG_INDEX j=0;j<R;j++){
			INDEX idx_1,idx_2;
			idx_1=j(11,12 - L_BRAMNUM);
			idx_2=j(12 - L_BRAMNUM - 1,0);
			out[j*C+position]=res[idx_1][idx_2];
		}

}

void step_4_para(in_type* out,long_uint res[2][bramnum][bramsize],uint32_t position){
		//transpose the matrix
		//write to position-th and (position+1-th) column
		typedef ap_uint<12> INDEX;
		typedef ap_uint<32> BIG_INDEX;
		step_4:for(BIG_INDEX j=0;j<R;j++){
			INDEX idx_1,idx_2;
			idx_1=j(11,12 - L_BRAMNUM);
			idx_2=j(12 - L_BRAMNUM - 1,0);
			in_type temp=res[0][idx_1][idx_2];
			temp=temp<<BASE_WIDTH;
			temp+=res[1][idx_1][idx_2];
			out[(j*C+position)/2]=temp;
		}

}

void Four_step_NTT_2_24_Ncore_large_bw(in_type* in,in_type* out,in_type* temp1,in_type* temp2){
	

	long_uint w_i = ONE;
	long_uint w_1=OMEGA[24];
	long_uint pre_w[4096];

	typedef ap_uint<12> INDEX;
	typedef ap_uint<32> BIG_INDEX;

	long_uint buffer_1[bramnum][bramsize];
#pragma HLS BIND_STORAGE variable=buffer_1 type=ram_t2p impl=uram
#pragma HLS ARRAY_PARTITION variable = buffer_1 complete dim = 1

//second copy
    long_uint buffer_2[bramnum][bramsize];
#pragma HLS BIND_STORAGE variable=buffer_2 type=ram_t2p impl=uram
#pragma HLS ARRAY_PARTITION variable = buffer_2 complete dim = 1

        //third copy
    long_uint buffer_3[bramnum][bramsize];
#pragma HLS BIND_STORAGE variable=buffer_3 type=ram_t2p impl=uram
#pragma HLS ARRAY_PARTITION variable = buffer_3 complete dim = 1


	//precompute w
	precompute:for(int i=0;i<C;i++){
#pragma HLS PIPELINE
		pre_w[i]=w_i;
		w_i = MUL_MOD(w_i, w_1,p);
	}

	//step-1
	step_1and2:for(int i=0;i<R;i+=3){
#pragma HLS UNROLL

		uint32_t position=i;
		//read i-th column
		read_column_reverse(in,buffer_1,position);

		// NTT_2_12_in_place(step1_in);
		ntt_2_12_2core_large_bw(buffer_1);

		//step-2
		step_2(temp1,buffer_1,pre_w[position],position);
		// w_i = MUL_MOD(w_i, w_1,p);

//second copy
		position+=1;
		if(position>=R)
            break;
		//read i-th column
		read_column_reverse(in,buffer_2,position);

		// NTT_2_12_in_place(step1_in);
		ntt_2_12_2core_large_bw(buffer_2);

		//step-2
		step_2(temp1,buffer_2,pre_w[position],position);
		// w_i = MUL_MOD(w_i, w_1,p);

		//third copy
		position+=1;
		if(position>=R)
            break;
		//read i-th column
		read_column_reverse(in,buffer_3,position);

		// NTT_2_12_in_place(step1_in);
		ntt_2_12_2core_large_bw(buffer_3);

		//step-2
		step_2(temp1,buffer_3,pre_w[position],position);
		// w_i = MUL_MOD(w_i, w_1,p);
	}

	step_3and4:for(int i=0;i<C;i+=3){
		#pragma HLS UNROLL
		uint32_t position=i;

		read_row_reverse(temp1,buffer_1,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_2core_large_bw(buffer_1);
		//step-4
		step_4(out,buffer_1,position);

		position+=1;
		if(position>=C)
            break;
		//second copy
		read_row_reverse(temp1,buffer_2,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_2core_large_bw(buffer_2);
		//step-4
		step_4(out,buffer_2,position);

		position+=1;
		if(position>=C)
            break;
		//second copy
		read_row_reverse(temp1,buffer_3,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_2core_large_bw(buffer_3);
		//step-4
		step_4(out,buffer_3,position);

	}
}

void Test(in_type* in,in_type* out){
	#pragma HLS BIND_STORAGE variable=OMEGA type=rom_np impl=lutram
	long_uint w_i = ONE;
	long_uint w_1=OMEGA[24];
	long_uint w_2=OMEGA[23];
	// long_uint pre_w[4096];

	typedef ap_uint<12> INDEX;
	typedef ap_uint<32> BIG_INDEX;

	long_uint step1_in[2][bramnum][bramsize];
//#pragma HLS BIND_STORAGE variable=step1_in type=ram_t2p impl=bram
#pragma HLS ARRAY_PARTITION variable = step1_in complete dim = 1
#pragma HLS ARRAY_PARTITION variable = step1_in complete dim = 2

//second copy
    long_uint step1_in_2[2][bramnum][bramsize];
//#pragma HLS BIND_STORAGE variable=step1_in_2 type=ram_t2p impl=bram
#pragma HLS ARRAY_PARTITION variable = step1_in_2 complete dim = 1
#pragma HLS ARRAY_PARTITION variable = step1_in_2 complete dim = 2

        //third copy
    long_uint step1_in_3[2][bramnum][bramsize];
//#pragma HLS BIND_STORAGE variable=step1_in_3 type=ram_t2p impl=bram
#pragma HLS ARRAY_PARTITION variable = step1_in_3 complete dim = 1
#pragma HLS ARRAY_PARTITION variable = step1_in_3 complete dim = 2

	//input size:R rows,C columns
	
	//step-1
	//int count=R;
    int count=4096;
	step_1and2:for(int i=0;i<count;i+=6){
#pragma HLS UNROLL

		uint32_t position=i;
		//read i-th column
		read_column_reverse_para2(in,step1_in,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_4core_large_bw(step1_in[0]);
        ntt_2_12_4core_large_bw(step1_in[1]);
		//step-2
		step_2_para(out,step1_in,w_i,position);
        w_i = MUL_MOD(w_i, w_2,p);
		// w_i = MUL_MOD(w_i, w_1,p);


		position+=2;
        if(position>=count)
            break;
		//read i-th column
		read_column_reverse_para2(in,step1_in_2,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_4core_large_bw(step1_in_2[0]);
        ntt_2_12_4core_large_bw(step1_in_2[1]);
		//step-2
		step_2_para(out,step1_in_2,w_i,position);
        w_i = MUL_MOD(w_i, w_2,p);
		// w_i = MUL_MOD(w_i, w_1,p);


		position+=2;
        if(position>=count)
            break;
		//read i-th column
		read_column_reverse_para2(in,step1_in_3,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_4core_large_bw(step1_in_3[0]);
        ntt_2_12_4core_large_bw(step1_in_3[1]);
		//step-2
		step_2_para(out,step1_in_3,w_i,position);
        w_i = MUL_MOD(w_i, w_2,p);
		// w_i = MUL_MOD(w_i, w_1,p);


	}

		count=C;

		step_3and4:for(int i=0;i<count;i+=6){
#pragma HLS UNROLL

		uint32_t position=i;
		//read i-th row
		read_row_reverse_para2(in,step1_in,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_4core_large_bw(step1_in[0]);
        ntt_2_12_4core_large_bw(step1_in[1]);
		//step-2
		step_4_para(out,step1_in,position);
		// w_i = MUL_MOD(w_i, w_1,p);


		position+=2;
        if(position>=count)
            break;
		//read i-th row
		read_row_reverse_para2(in,step1_in_2,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_4core_large_bw(step1_in_2[0]);
        ntt_2_12_4core_large_bw(step1_in_2[1]);

		step_4_para(out,step1_in_2,position);


		position+=2;
        if(position>=count)
            break;
		//read i-th row
		read_row_reverse_para2(in,step1_in_3,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_4core_large_bw(step1_in_3[0]);
        ntt_2_12_4core_large_bw(step1_in_3[1]);

		step_4_para(out,step1_in_3,position);

	}
}


void Test_Uram(in_type* in,in_type* out){
	long_uint w_i = ONE;
	long_uint w_1=OMEGA[24];
	long_uint w_2=OMEGA[23];
	// long_uint pre_w[4096];

	typedef ap_uint<12> INDEX;
	typedef ap_uint<32> BIG_INDEX;

	long_uint buffer_0[2][bramnum][bramsize];
#pragma HLS BIND_STORAGE variable=buffer_0 type=ram_t2p impl=uram
#pragma HLS ARRAY_PARTITION variable = buffer_0 complete dim = 1
#pragma HLS ARRAY_PARTITION variable = buffer_0 complete dim = 2

//second copy
    long_uint buffer_1[2][bramnum][bramsize];
#pragma HLS BIND_STORAGE variable=buffer_1 type=ram_t2p impl=uram
#pragma HLS ARRAY_PARTITION variable = buffer_1 complete dim = 1
#pragma HLS ARRAY_PARTITION variable = buffer_1 complete dim = 2

        //third copy
    long_uint buffer_2[2][bramnum][bramsize];
#pragma HLS BIND_STORAGE variable=buffer_2 type=ram_t2p impl=uram
#pragma HLS ARRAY_PARTITION variable = buffer_2 complete dim = 1
#pragma HLS ARRAY_PARTITION variable = buffer_2 complete dim = 2

	//input size:R rows,C columns
	
	//step-1
	//int count=R;
    int count=64;
	step_1and2:for(int i=0;i<count;i+=6){
#pragma HLS UNROLL

		uint32_t position=i;
		//read i-th column
		read_column_reverse_para2(in,buffer_0,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_4core_large_bw(buffer_0[0]);
        ntt_2_12_4core_large_bw(buffer_0[1]);
		//step-2
		step_2_para(out,buffer_0,w_i,position);
        w_i = MUL_MOD(w_i, w_2,p);
		// w_i = MUL_MOD(w_i, w_1,p);


		position+=2;
        if(position>=count)
            break;
		//read i-th column
		read_column_reverse_para2(in,buffer_1,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_4core_large_bw(buffer_1[0]);
        ntt_2_12_4core_large_bw(buffer_1[1]);
		//step-2
		step_2_para(out,buffer_1,w_i,position);
        w_i = MUL_MOD(w_i, w_2,p);
		// w_i = MUL_MOD(w_i, w_1,p);


		position+=2;
        if(position>=count)
            break;
		//read i-th column
		read_column_reverse_para2(in,buffer_2,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_4core_large_bw(buffer_2[0]);
        ntt_2_12_4core_large_bw(buffer_2[1]);
		//step-2
		step_2_para(out,buffer_2,w_i,position);
        w_i = MUL_MOD(w_i, w_2,p);
		// w_i = MUL_MOD(w_i, w_1,p);


	}

		//count=C;

		step_3and4:for(int i=0;i<count;i+=6){
#pragma HLS UNROLL

		uint32_t position=i;
		//read i-th row
		read_row_reverse_para2(in,buffer_0,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_4core_large_bw(buffer_0[0]);
        ntt_2_12_4core_large_bw(buffer_0[1]);
		//step-2
		step_4_para(out,buffer_0,position);
		// w_i = MUL_MOD(w_i, w_1,p);


		position+=2;
        if(position>=count)
            break;
		//read i-th row
		read_row_reverse_para2(in,buffer_1,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_4core_large_bw(buffer_1[0]);
        ntt_2_12_4core_large_bw(buffer_1[1]);

		step_4_para(out,buffer_1,position);



		position+=2;
        if(position>=count)
            break;
		//read i-th row
		read_row_reverse_para2(in,buffer_2,position);
		// NTT_2_12_in_place(step1_in);
		ntt_2_12_4core_large_bw(buffer_2[0]);
        ntt_2_12_4core_large_bw(buffer_2[1]);

		step_4_para(out,buffer_2,position);


	}
}

void test_para(in_type* in1,in_type* out1,in_type* in2,in_type* out2){
    Test(in1,out1);
    Test(in2,out2);
	// Test(in3,out3);
    // Test(in4,out4);
//    Test_Uram(in5,out5);
//    Test_Uram(in6,out6);

}
