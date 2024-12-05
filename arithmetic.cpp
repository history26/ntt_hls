#include "arithmetic.hpp"

using namespace std;

long_uint add_long(long_uint left,long_uint right){
	return left-right;
}
long_uint sub_long(long_uint left,long_uint right){
	return left-right;
}

//res=left*right
long_long_uint mult_long(long_long_uint left,long_long_uint right){
	long_long_uint res=left*right;
	return res;
}

//res = (left*right) mod b
long_uint mult_mod_base_long(long_long_uint left,long_long_uint right){

	long_long_uint res=mult_long(left,right);
	long_uint res_width=res(BASE_WIDTH-1,0);
	return res_width;
}

long_uint mul_partly(long_uint left,long_uint right){
	long_uint ll=left(127,0);
	long_uint lh=left(255,128);
	long_uint rl=right(127,0);
	long_uint rh=right(255,128);

	long_uint mid=(ll*rh+lh*rl)<<128;
	long_uint res=mid+ll*rl;
	return res;

}

long_long_uint mul_long_opt(long_uint left,long_uint right){
	long_uint ll=left(127,0);
	long_uint lh=left(255,128);
	long_uint rl=right(127,0);
	long_uint rh=right(255,128);

	long_uint HH=lh*rh;
	long_uint LL=ll*rl;
	long_uint mid=(lh+ll)*(rh+rl)-HH-LL;

	long_long_uint res=(long_long_uint(HH)<<BASE_WIDTH)+(long_long_uint(mid)<<128)+LL;
	return res;

}


//res = (-x) mod m
long_uint nevigate(long_uint x,long_uint m){
	long_uint res =sub_long(m,x);
	return res;
}


//return t*R^(-1) mod m
//R = b^n,b = 2^LONG_WIDTH,n = 1,negat_m_inv = -m^(-1) mod b
long_uint montgomery_reduce(long_long_uint t,long_uint m,long_uint INV){
	//m'=-m^(-1) mod b
	long_uint t_mod_b=t>>BASE_WIDTH;
	long_uint u=mul_partly(t_mod_b,INV);
	long_long_uint temp=mul_long_opt(u,m);
	long_long_uint a=t+temp;
	long_uint res=a>>BASE_WIDTH;
	if(res>=m){
		res=res-m;
	}
	return res;
}

//res = left * right ,then montgomery_reduce
long_uint MUL_MOD(long_long_uint left,long_long_uint right,long_uint m){
	long_long_uint temp=mul_long_opt(left,right);
	long_uint res=montgomery_reduce(temp,m,INV);
	return res;
}

//res = (left - right) mod m , underflow will be resolve, left and right is smaller than m
long_uint SUB_MOD(long_uint left,long_uint right,long_uint m){
	long_uint tmp=m-right;
	long_uint res=left+tmp;
	if(res>=m){
		res=res-m;
	}
	return res;
}

//res = left + right , makesure overflow will be resolve
long_uint ADD_MOD(long_uint left,long_uint right,long_uint m){
	long_uint res=left+right;
	if(res>=m){
		res=res-m;
	}
	return res;
}
