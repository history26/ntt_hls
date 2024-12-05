#define corenum 4
#define R 4096
#define C 4096
// #define N 16777216
#define N 4096
#define IO_para 2
#define stagemax 10
#define STAGENUM 12

#if (corenum==1)
	#define bramnum 2
	#define bramsize 2048
	#define L_BRAMNUM 1
	#define L_BRAMSIZE 11
#elif (corenum==2)
	#define bramnum 4
	#define bramsize 1024
	#define L_BRAMNUM 2
	#define L_BRAMSIZE 10
#elif (corenum==4)
	#define bramnum 8
	#define bramsize 512
	#define L_BRAMNUM 3
	#define L_BRAMSIZE 9
	#define RPBRAMNUM 4
	#define RPBRAMSIZE 1024
	#define L_RPBRAMNUM 2
	#define L_RPBRAMSIZE 10
#elif (corenum==8)
	#define bramnum 16
	#define bramsize 256
	#define L_BRAMNUM 4
	#define L_BRAMSIZE 8
#elif (corenum==16)
	#define bramnum 32
	#define bramsize 128
	#define L_BRAMNUM 5
	#define L_BRAMSIZE 7
#endif
