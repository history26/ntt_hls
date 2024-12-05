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

#include "top.h"
#include <iostream>
#include <fstream>

using namespace std;

//result is large unsigned
//high->low:in[0],in[1],鈥︹��,in[size-1]
long_uint u64_2_ularge(uint64_t *in,int size){
	long_uint tmp=0;
	tmp+=in[0];
	for(int i=1;i<size;i++){
		tmp=tmp<<64;
		tmp+=in[i];
	}
	return tmp;
	
}

//result is large signed
//high->low:in[0],in[1],鈥︹��,in[size-1]
long_int u64_2_large(uint64_t *in,int size){
	long_int tmp=(long_int)in[0];
	for(int i=1;i<size;i++){
		tmp=(tmp<<64)+(long_int)in[i];
	}
	return tmp;
	
}


int ntt_test(){

	//example input and output
	static in_type test_in[R/2];
	static in_type test_out[R/2];
	// long_uint* test_in=(long_uint*)malloc(N);
	// long_uint* test_out=(long_uint*)malloc(N);

	static in_type temp1[R/2];
	static in_type temp2[R/2];

	static long_uint true_out[R];

	string datapath="./data/";
    string infile=datapath+"input.txt";
    string outfile=datapath+"output.txt";
	int err=0;
	int size=4;

    fstream instream(infile.c_str(),ios::in);
    fstream outstream(outfile.c_str(),ios::in);

	if(instream.is_open()&&outstream.is_open()){

    	for(uint64_t i=0;i<R;i++){
			uint64_t t1[size];
			uint64_t t2[size];
			for(int j=size-1;j>-1;j--){
//				fscanf("%llx",t1[j]);
				instream>>t1[j];
				outstream>>t2[j];
			}
			test_in[i]=u64_2_ularge(t1,size);

			true_out[i]=u64_2_ularge(t2,size);
    	}
    }
    else{
    	cout<<"in/out file is not open"<<endl;
    }
   instream.close();
   outstream.close();



	// ap_uint<64> test[4][8]={{0x0D015766AFF0868B,0x0D015766AFF0868B,0x0D015766AFF0868B,0xBBBB67F06FF5E5EC},
	// 	{0xFC1988E68A9F75E0,0xFC1988E68A9F75E0,0xFC1988E68A9F75E0,0x3CF52321D4F161D5},
	// 	{0x0ECBE6CA5AF3868E,0x0ECBE6CA5AF3868E,0x0ECBE6CA5AF3868E,0xDD3775487E5B9104},
	// 	{0x0C1BA4B4EB5DC119,0x0C1BA4B4EB5DC119,0x0C1BA4B4EB5DC119,0xCC4AA6CA52FE1687}};

	// for(int i=0;i<N;i++){
	// 	test_in[i]=ZERO;
	// }
	// test_in[0]=ONE;
	// test_in[R]=TWO;
	dut(test_in,test_out,temp1,temp2);

	// GF expected_out[8] = { GF(3ull), GF(16539086126023702892ull), GF(
	// 	911812899281015199ull), GF(2932311101274590780ull), GF(4259936512344175332ull),
	// 	GF(3503108039395815018ull), GF(15858278158945469571ull), GF(
	// 			13920042588279902734ull) };


	for(int i=0;i<8;i++){

		cout << "NTT result [" << i << "] should be "
				<< test_out[i].to_string(16) << ", but the output is "
				<< true_out[i].to_string(16) << endl;

	}
	
	return 0;

}


int main() {

	int res=ntt_test();
	
	
	return 0;
}
