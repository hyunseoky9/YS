#define _block_ 30
#define _PRINTOUT_EVERY_ 1000
#define _DELTA_ 0.1
#define _DELTAover2_ 0.05
#define _haplorank_ 5

#define noStorageMode 0
#define multilocus 1



//print options for visualizing : on - off option
#define PRINT_WHERE 0
#define PRINT_WHERE_SHORT 0
#define PRINT_RESULT 0
#define PRINT_TML 0				//PRINT_ASSIGN_NONSYN & PRINT_SELOPT & PRINT_DELTA
#define PRINT_RECOMB 0
#define PRINT_SELECTION 0
#define PRINT_GETFREQ 0
#define PRINT_FIT 0
#define PRINT_SEQ 0

//on-off option 
#define TRJ 1
#define DELTA_FIXED 0
#define _SearchingPara_ 0





/*========================================IMGSOB Version Info==============================================


	//printf("\n%d", __LINE__);	getchar(); 
	//useful tool for debugging
	
	//getFrequency 바꾸기 -> 각 자리에 대해 만 제너레이션까지 더하고, SE구할 수 있게 만들기
	
		//매 만 제너레이션마다 멤셋해주고
		//이닛시뮬레이션에서 멜럭해주기(엔시퀀스길이만큼)
		//포인터로 리턴받기
		
	//한 제너레이션에서 다음 제너레이션으로 잘 넘어가고 있나 보기
	
	//셀렉션옵티마 모드를 인풋파일에 넣을까?	
	
	//시그마 설정 어떻게?? 지금은 그냥 1임
	
	//임시변수선언해둔것들 몽땅 해결하기
	
	//시드는 전역 받아서 뽑는걸로 다 통일하기!!
	
	// 모두 완료하면 neutrality test

	
	//0330		focused on <<seachingParameter>> which enables balancing selection		
					there are only one nonsyn site, one mutant individual as default
					check lossOrfix to find balancing selection trial
					get the longest generation at the end of the simulation
					if balancing selection(gen = nGen) does not occur until (nNonCountTry==10000), printf("No balancing selection");
					no recombination
					
	//0331			_DELTAover2_ defined : two fitness graph should be _DELTAover2_를 중심으로 symmetrical(0330 notes)
					
	
	//0404		parameters rescaling (0403 notes)	
					- miminum과 maximum selOpt은 0을 중심으로 Delta(큰델타. delta아님)/2 씩 떨어져 있다
						-> _DELTAover2_ no longer used
					- 매 trial마다 pheno_wild는 -sMax와 sMax사이의 uniform distribution에서 랜덤하게 설정
						+ 더 이상 DELTA_FIXED 옵션을 사용하지 않음
					- delta ~ Normal(0, sMax)
				
				No more searchingParameter
					
				inputFile 	->	20 nonsynSites
	
	//0405		GetHet	- using GetHet instead of GetFrequency
					- FreqMut is divided by (nIndT)^2 *_PRINTOUT_EVERY_
	
	
	//0407		Nmut 잘못되어있었음. p.nSeq_block으로 바꿨음
	
	//0411		0410 노트 참조
					해당 generation에서 가장 빈번한 haplotype을 구하는 haplotypeTable을 만들기 위해 
					연결리스트 구현
					
					-> SvHPF.c
				
	//0424		svardal_afterdebug_Guli
	
	//0425		svardal_afterdebug_Guli_sameFit
					//
					
	----------------------------svardalGul_0515.c-----------------------------------------------------------	
	//0515		- del_range가 parameter 에 추가됨	=> delta ~ Normal (0, del_range)
	//0516		- print_Haplo_outfile : top5개만 print하는 것이 아니라 전체 haplotype 프린트
					gen0에서도 프린트		=> if ((gen % _PRINTOUT_EVERY_ == _PRINTOUT_EVERY_-1) || (gen== 0))
				
				- print delta in out file
				
				- 하나 뮤테이션 자동으로 주고 있었음 (g_rSite)
					-> 이것이 그래프의 패턴이나 invasion과 관련이 있을까?
			
				
	----------------------------svardalGul_0531.c-----------------------------------------------------------	
	<SeedBank Model> 
	 - new fuctions : Mutation_SeedBank, SelectionAndReproducation_Guli_SeedBank
		- Mutation_SeedBank 					  
			: no mutation in second subpop
		- SelectionAndReproducation_Guli_SeedBank 
			: in second population, there is no genetic drift but the order of individual is randomly sorted
			
	//0602	free memories, especially haploTable and ranArray
	//0619	print_Haplo_outfile in generation 0 is fixed -> before this version, haplo type is printed every generation
	
	//0904	delta의 분포를 gasdev 에서 expdev로 고침. (AssignDeltaEachNon_Exp)
	
	----------------------------svardalGul_0908_7.c------------------------------------------------------------
	<SeedBank + Recombination Model>
	 - 0905 ~ 0908_6 : don't use. 0908_7.c is the debugged version.
	 - modified functions : 
	 
	 
	//0913	delta의 분포를 다시 expdev에서 gasdev로 고침 (AssignDeltaEachNon_Exp)
	
	----------------------------svardalGul_0918.c------------------------------------------------------------
	<.freq + .haplo + .out>
	 - 더 이상 .pheno를 뽑지 않음 -> .haplo를 every generation에서 뽑는 것으로 대체
	 * modified fuctions 	
		: main	-	[v] haploCount의 위치가 main()에서 SelectionAndReproducation_Guli_SeedBank()로 이동
					[v] 초기화한 population에서의 haplotype은 뽑지 않는다 (이전 버전의 gen0에 해당)
					[v] phenotype에 인쇄하는 부분 모두 주석처이
					
		: SelectionAndReproducation_Guli_SeedBank()
				- [v] phenotype인쇄부분 삭제
				- [v] recombination 이후 (CopyIndividual바로 뒤)에 haploCount를 넣음. 
		
		: Recombination()	-	[v] recomRate == 0일 경우 바로 return ; no recom simulation을 위해
		
	 
	 * 에러 수정 (진행상황 : 코드만 고침 - 결과 확인 필요 )
		- main에서 if ( periodCount == 0) -> if ( periodCount%period == 0)
		 : if ( periodCount == 0)로 하면 DisturbSelOpt이 gen = 0에서만 실행된다. 
			이렇게 진행되어왔다면 사실상 disturbSelOpt이 작동하지 않은 거나 다름없는 수준. 
			실제 결과에서 Distrub옵션이 작동했을 때와 아닐 때의 결과가 비슷한지 확인해야 한다. 
	
	
	----------------------------svardalGul_1016.c------------------------------------------------------------
	<back to refuge model>
	- 컴파일시 SG_DRDR_SEEDR말고 SG_DFDR_REFUR로 컴파일할 것
	 * new functions
		: SelectionAndReproducation_Guli_Refuge
		
	 * modified functions
		: main	-	[v] mutation_SeedBank를 mutation으로 교체
				-	[v] SelectionAndReproducation_Guli_SeedBank 를 SelectionAndReproducation_Guli_Refuge 로 교체
		
		: mutation
				-	[v] Nu * 2 해주기 -> mutation_SeedBank에서는 subpop1에만 mutation이 일어남 
											mutation에서는 1, 2에서 모두 일어나야 하므로 N이 두배임
											우리는 Nu 자체를 파라미터로 받고 있으므로 수동으로 *2 해줘야함
============================================================================================================*/



#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>			//standardized bit 
#include <string.h>
#include <omp.h>
#include <math.h>
#define PI 3.141592654

//  
#if defined _SBE_
    #define _MUT1_
	#define _SELSB_
	
#elif defined _SBC_
	#define _MUT12_
    #define _SELC_
	
#elif defined _RSE_
	#define _MUT12_
	#define _SELRS_
	
#elif defined _RSC_
	#define _MUT12_
	#define _SELC_
	

#endif  

//global variables
FILE *g_outfile, *g_infile, *g_freqfile, *g_phenofile, *g_haplofile;
long gseed;
int g_rSite;


struct Individual{
	unsigned int* seq;
	
};

struct Population {
	struct Individual* indiv;
};

typedef struct {
	int nDeme;
	int nTry;
	int nGen;
	
	int nInd1;
	int nInd2;
	int nIndT;
	int* popsize;
	
	int nSeq_block;
	int nSeq;
	int nNon;
	int nSyn;
	
	int mig;
	int period;
	
	double Nu;
	double recomRate;
	
	double sMax;
	double epsilon;
	double C;
	double sigma_s;
	double del_range;
	
	long seed;
	
}para_t, *Ptrpara_t;


//linked list for ranking haplotype - modified linked list including int array for sequece info
typedef struct _node
{
	struct _node *next;
	unsigned int* seq;
	int numSeq;
}node;


typedef node* nptr;

typedef struct _list
{
	int count;
	nptr head;
}list;


//Printing fuctions------------------------------------
//void PrintBinary(unsigned int input);
//void PrintBinary_out(unsigned int input, FILE* file );
//void PrintSeq_non_out(unsigned int* seq,unsigned int* nonsyn, FILE* file , para_t p);
//void PrintBinary_asterik(unsigned int input);
//void PrintPop(struct Population* pop, para_t p_);

//void init(list* lptr);
//void insert(list* lptr, int numSeq, unsigned int* seq, int position, para_t p);
//void print_top5Haplo_outfile(list* lptr, unsigned int* nonsyn, double* delta, double pheno_wild, para_t p);
//void countHaplotype(list*lptr, unsigned int* seq, para_t p)			//countHaplotype


//random fuctions---------------------------------
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI 3.141592654

//float ran1(long *idum);
//float gammln(float xx);
//float poidev(float xm, long *idum);
//float bnldev(float pp, int n, long *idum);
//float gasdev(long *idum);


//initializing functions----------------------
//para_t ReadParameter(const char* inpname, const char* outname,  para_t p);
//void AssignNonsynSite(unsigned int* template_nonsys, para_t p);
//void AssignDeltaEachNon(double* delta, para_t p);
//void AssignSelOptima(double** selOpt, int mode, para_t p);
//struct Population* allocPop(para_t p);


//struct Individual CopyIndividual(struct Individual dest, struct Individual sourc, para_t p);

//double getPhenoFromSeq(unsigned int* seq, unsigned int* nonsyn, double* delta, double pheno_wild, para_t p);
//double FitOfIndiv(struct Individual* ind, unsigned int* nonsyn, double* delta, double pheno_wild, double selOpt, para_t p);

//void Mutation(struct Population* pop, para_t p);
//void Migration(struct Population* pop, para_t p);
//void Recombination(struct Individual* ind1, struct Individual* ind2, para_t p);
//struct Population* SelectionAndReproducation(struct Population* currentPop, struct Population* pop1, struct Population* pop2,unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen, para_t p);
//void PrintFreq(struct Population* pop, int gen, para_t p);
//int PrintFreqAndFixLossCheck(struct Population* pop, int gen, para_t p);

void PrintBinary(unsigned int input)	
{
	int i, check = 1;
	char bit;
	
	for (i=_block_-1; i>=0; i--)
	{
		bit = (input & (1 << i))?1:0;
		//if (check)
		//{
			if (bit ==1)
			{
				check = 0;
				printf("%d", bit);
			}
		//}
		else printf("%d", bit);
	}
	//printf("\n");
}

void PrintBinary_out(unsigned int input, FILE* file )	//has changed
{
	int i, check = 1;
	char bit;
	
	for (i=_block_-1; i>=0; i--)
	{
		bit = (input & (1 << i))?1:0;
		//if (check)
		//{
			if (bit ==1)
			{
				check = 0;
				fprintf(file,"%d", bit);
			}
		//}
		else fprintf(file,"%d", bit);
	}
	//printf("\n");
}


void PrintSeq_non_out(unsigned int* seq,unsigned int* nonsyn, FILE* file , para_t p)	//has changed
{
	int i, j,  check = 1;
	char bit;
	unsigned int seg, template_seg;
	int count = 0;
	
	
	for (j = 0; j < p.nSeq_block; j++)
	{
		seg = seq[j];
		template_seg = nonsyn[j];
		count = 0;
		
		for (i=_block_-1; i>=0; i--)
		{
			
			
			if ((template_seg & (1 << i))?1:0)			//if this site is nonsyn site
			{
				//printf("block %d %d %d %lf", i,  count, delta_idx, z);
				bit =(seg & (1 << i))? 1: 0;
	
				if (bit ==1)
				{
					check = 0;
					fprintf(file,"%d", bit);
				}
				else 
					fprintf(file,"%d", bit);
				
				
		
			}
			
		}
		
		fprintf(file,"/");
	}
	
				
	
	//printf("\n");
}

void PrintSeq_all_out(unsigned int* seq,unsigned int* nonsyn, FILE* file , para_t p)	//has changed
{
	int i, j,  check = 1;
	char bit;
	unsigned int seg, template_seg;
	int count = 0;
	
	
	for (j = 0; j < p.nSeq_block; j++)
	{
		seg = seq[j];
		template_seg = nonsyn[j];
		count = 0;
		
		for (i=_block_-1; i>=0; i--)
		{
			
			
			
				//printf("block %d %d %d %lf", i,  count, delta_idx, z);
				bit =(seg & (1 << i))? 1: 0;
	
				if (bit ==1)
				{
					check = 0;
					fprintf(file,"%d", bit);
				}
				else 
					fprintf(file,"%d", bit);
				
				
		
			
			
		}
		
		fprintf(file," ");
	}
	
				
	
	//printf("\n");
}




void PrintBinary_asterik(unsigned int input)	
{
	int i, check = 1;
	char bit;
	
	for (i=_block_-1; i>=0; i--)
	{
		bit = (input & (1 << i))?1:0;
		//if (check)
		//{
			if (bit ==1)
			{
				check = 0;
				printf("*");
			}
		//}
		else printf("-");
	}
	//printf("\n");
}

void PrintPop(struct Population* pop, para_t p_)
{
	int dm, i, s;
	//int count = 0;
	int seg_print;
	
	
	for (dm = 0; dm<p_.nDeme; dm++)
	{
		
		printf(" > SubPopulation %d \n", dm);
		for (i = 0 ; i < p_.popsize[dm] ; i++)
		{
			//sIdx_check = 0;	
			printf("indiv %d : ",  i);
			
			for (s = 0; s < p_.nSeq_block; s++)
			{
				if(s==0)		
					printf("%p ", &(pop[dm].indiv[i].seq[0]));
				
				seg_print = pop[dm].indiv[i].seq[s];
				PrintBinary(seg_print);
				//printf("(%d)", seg_print);
				printf(" ");
				//count++;
				
			}
			printf("\n");
		}
	}
	printf("\n\n");
	
}


void init(list* lptr)
{
	//initialize the list
	if (lptr->head != NULL)
	{
		nptr tmp1 = lptr->head;
		nptr tmp2;
		
		while (tmp1 != NULL)
		{
			tmp2 = tmp1;
			tmp1 = tmp2->next;
			free(tmp2->seq);
			free(tmp2);
			
		}
			
	}
	
	lptr->count = 0;
	lptr->head = NULL;
	
	
		
	
}

void insert(list* lptr, int numSeq, unsigned int* seq, int position, para_t p)
{

	//insert value to the proper postion
	if ( position < 1 || position > (lptr->count)+1 )
	{
		printf("Position Out of Bound\n");
		return;
	}
	
	
	nptr new_nptr = (node*)malloc(sizeof(node));
	new_nptr->numSeq = numSeq;
	new_nptr->seq = (unsigned int*)malloc(sizeof(unsigned int)*p.nSeq_block);
	
	memcpy(new_nptr->seq, seq, sizeof(unsigned int)*p.nSeq_block);
	
	/*int i; //
	
	printf("> ");
	for (i=0; i<p.nSeq_block; i++)
	{
		printf("%d ", new_nptr->seq[i]);
	}
	printf("\n>>");
	for (i=0; i<p.nSeq_block; i++)
	{
		printf("%d ", seq[i]);
	}
	printf("\n");
	*/
	
	if (position == 1)
	{
		new_nptr->next = lptr->head;
		lptr->head = new_nptr;
	}
	else
	{
		nptr tmp = lptr->head;
		int i;
		for (i = 1; i <position-1; i++)
		{
			tmp = tmp->next;
		}
		new_nptr->next = tmp->next;
		tmp->next = new_nptr;
	}
	lptr->count++;
}

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}


double getPhenoFromSeq(unsigned int* seq, unsigned int* nonsyn, double* delta, double pheno_wild, para_t p)
{
	
	int i, count, delta_idx=0;
	double z = 0;
	unsigned int seg, template_seg;
											
	for (i = 0; i < p.nSeq_block; i++)
	{
		count = 0;
		seg = seq[i];
		template_seg = nonsyn[i];
		
		
		while ( count < _block_)
		{
			if ((template_seg & (1 << count))?1:0)			//if this site is nonsyn site
			{
				//printf("block %d %d %d %lf", i,  count, delta_idx, z);
				if ((seg & (1 << count))? 1: 0)				//if the allele is mutant type
				{	
					//printf(" V");
					z+= delta[delta_idx];
					
				}
				//printf("\n");
				delta_idx++;
					
			}
			count++;
			
		}
		
		
	}
	//printf("z : %lf\n", z);
	
	z += pheno_wild;
	
	return z;
	
	
}


void print_top5Haplo_outfile(list* lptr, unsigned int* nonsyn, double* delta, double pheno_wild, para_t p)
{
	unsigned int* seq;
	double pheno;
	int i, rank=0;
	
	
	nptr tmp = lptr->head;
	while(tmp != NULL)
	//while(tmp->next != NULL)		//0516 corrected
	{
		seq = tmp->seq;
		 
		/* print all sites
		for (i=0; i<p.nSeq_block; i++)
		{
			PrintBinary_out(seq[i], g_haplofile);
			fprintf(g_haplofile, "/");
		}
		*/
		//fprintf(g_haplofile, ">");
		//PrintSeq_non_out( seq, nonsyn, g_haplofile , p);
		PrintSeq_all_out( seq, nonsyn, g_haplofile , p);
		//fprintf(g_haplofile, "<");
		
		/*
		for (i=0; i<p.nSeq_block; i++)
		{
			fprintf(g_haplofile, "%d ", seq[i]);
		}
		*/
		pheno = getPhenoFromSeq(seq, nonsyn, delta, pheno_wild, p);
		 
		
		fprintf(g_haplofile, "%lf ", pheno );
		 
		
		fprintf(g_haplofile, "%d\n", tmp->numSeq);
		 
		
		
		rank++;
		if (rank > _haplorank_) 
			break;
		tmp = tmp->next;
		
		
	}
	
}


void print_Haplo_outfile(list* lptr, unsigned int* nonsyn, double* delta, double pheno_wild, para_t p)
{
	unsigned int* seq;
	double pheno;
	int i, rank=0;
	
	
	nptr tmp = lptr->head;
	while(tmp != NULL)
	//while(tmp->next != NULL)		//0516 corrected
	{
		seq = tmp->seq;
		 
		/* print all sites
		for (i=0; i<p.nSeq_block; i++)
		{
			PrintBinary_out(seq[i], g_haplofile);
			fprintf(g_haplofile, "/");
		}
		*/
		//fprintf(g_haplofile, ">");
		//PrintSeq_non_out( seq, nonsyn, g_haplofile , p);
		PrintSeq_all_out( seq, nonsyn, g_haplofile , p);
		 
		//fprintf(g_haplofile, "<");
		
		/*
		for (i=0; i<p.nSeq_block; i++)
		{
			fprintf(g_haplofile, "%d ", seq[i]);
		}
		*/
		pheno = getPhenoFromSeq(seq, nonsyn, delta, pheno_wild, p);
		 
		
		fprintf(g_haplofile, "/%lf ", pheno );
		 
		
		fprintf(g_haplofile, "%d\n", tmp->numSeq);
		 
		
		
		rank++;
		//if (rank > _haplorank_) 
		//	break;
		tmp = tmp->next;
		
		
	}
	
}

void countHaplotype(list*lptr, unsigned int* seq, para_t p)			//countHaplotype
{
	//modified search()
	//transverse the list and 
	//find the first position of the value (first form head)
	//if not exist, return 0
	
	nptr tmp = lptr->head;
	int i = 1;
	int j;
	int same = 1;
	while ( tmp != NULL )
	{
		same = 1;
		for ( j = 0; j < p.nSeq_block; j++)
		{
			if (seq[j] != tmp->seq[j])
			{
				same = 0;
			}
		}
		
		if (same == 1 )
		{
			tmp->numSeq++;			
			break;
		}
		
		i++;
		tmp = tmp->next;
	}
	
	if (i > lptr->count )
	{
		insert(lptr, 1, seq, i,  p);
	}

}




//gammlm
float gammln(float xx)
{
	double x, y, tmp, ser;
	static double cof[6] = {76.18009173,-86.50532033,24.01409822,-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;
	
	y = x = xx;
	tmp = x +5.5;
	tmp -= (x+ 0.5) * log(tmp)	;
	ser = 1.000000000190015;
	for (j=0;j<=5;j++)	ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
	
//poidev 
#define PI 3.141592654
float poidev(float xm, long *idum)
{
	float gammln(float xx);
	float ran1(long *idum);
	static float sq, alxm, g, oldm=(-1.0);
	float em, t, y;
	
	if (xm < 12.0)
	{
		if (xm != oldm)
		{
			oldm = xm;
			g = exp(-xm);
		}
		em = -1;
		t = 1.0;
		do{
			++em;
			t *= ran1(idum);
		}
		while (t>g);
	}
	else
	{
		if (xm != oldm)
		{
			oldm = xm;
			sq = sqrt(2.0*xm);
			alxm = log(xm);
			g = xm*alxm - gammln(xm+1.0);
		}
			do {
			do {
					y = tan(PI*ran1(idum));
					em = sq*y+xm;
			}while (em < 0.0);
			em = floor(em);
			t = 0.9*(0.1+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		}while (ran1(idum)>t);
	}
	return em;
}


float bnldev(float pp, int n, long *idum)
{
	int j;
	static int nold=(-1);
	float am,em,g,angle,p,bnl,sq,t,y;
	static float pold=(-1.0),pc,plog,pclog,en,oldg;
	
	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (ran1(idum) < p) bnl += 1.0;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran1(idum);
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran1(idum);
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
				-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran1(idum) > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
	
}

float gasdev(long *idum)
{
	static int iset = 0;
	static float gset;
	float fac, rsq, v1, v2;
	
	if (*idum < 0) iset=0; 
	if (iset == 0) 
	{
		 do 
		 {
			v1=2.0*ran1(idum)-1.0; // pick two uniform numbers in the square extending from -1 to +1 in each direction, 
			v2=2.0*ran1(idum)-1.0; 
			rsq=v1*v1+v2*v2; // see if they are in the unit circle,
		}
		while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again. 
		fac=sqrt(-2.0*log(rsq)/rsq); /* Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time. */
		gset=v1*fac;
		iset=1; // Set ag. 
		return (double) v2*fac; 
	} 
	else 
	{ // We have an extra deviate handy, 
		iset=0; // so unset the ag, 
		return (double) gset; // and return it. 
	} 
}
		
		

float expdev(long *idum)
{
	float ran1(long *idum);
	float dum;

	do
		dum=ran1(idum);
	while (dum == 0.0);
	return -log(dum);
}

		
para_t ReadParameter(const char* inpname, const char* outname,  para_t p)
{

	int i;

	g_outfile = fopen(outname, "w");
	if ( (g_infile = fopen(inpname, "r")) == NULL)
    {
        printf("Cannot open input file\n");
        exit(-1);           //* <stlib.h> should be included. close every file when there is an error 
    }
	//read parameters from outFile
	
	fscanf(g_infile, "%*s%d", &(p.nDeme));
	fscanf(g_infile, "%*s%d", &(p.nTry));
	fscanf(g_infile, "%*s%d", &(p.nGen));
	
	/*fscanf(g_infile, "%*s%d", &(p.nInd1));
	fscanf(g_infile, "%*s%d", &(p.nInd2));
	p.nIndT = p.nInd1 + p.nInd2;
	*/
	
	p.nIndT = 0;
	p.popsize = (int*) malloc (sizeof(int) * p.nDeme);
	for (i =0 ; i < p.nDeme; i++)
	{
		fscanf(g_infile, "%*s%d", &(p.popsize[i]));
		//vi//printf("> %d : %d ", i, p.popsize[i]);
		p.nIndT += p.popsize[i];
	}
	
	
	
	fscanf(g_infile, "%*s%d", &(p.nSeq_block));
	fscanf(g_infile, "%*s%d", &(p.nNon));
	p.nSeq = p.nSeq_block * _block_;
	p.nSyn = p.nSeq - p.nNon;
	

	
	fscanf(g_infile, "%*s%d", &(p.mig));
	fscanf(g_infile, "%*s%d", &(p.period));

	fscanf(g_infile, "%*s%lf", &(p.Nu));
	fscanf(g_infile, "%*s%lf", &(p.recomRate));
	

	
	fscanf(g_infile, "%*s%lf", &(p.sMax));
	fscanf(g_infile, "%*s%lf", &(p.epsilon));
	fscanf(g_infile, "%*s%lf", &(p.C));
	fscanf(g_infile, "%*s%lf", &(p.sigma_s));
	fscanf(g_infile, "%*s%lf", &(p.del_range));
	
	fscanf(g_infile, "%*s%ld", &(p.seed));
	
	
	
	//print parameters to outFile
	#if defined(_SBC_)
	printf( " << SEEDBANK - CONTROL    MODE >> \n");
	#elif defined(_SBE_)
	printf( " << SEEDBANK - EXPERIMNET MODE >> \n");
	#elif defined(_RSC_)
	printf( " << REFUGE   - CONTROL    MODE >> \n");
	#elif defined(_RSE_)
	printf( " << REFUGE   - EXPERIMNET MODE >> \n");
	#else 
	printf("++++++++++++++++++++++++++++++\n+[Warning] Please compile with option+\n++++++++++++++++++++++++++++++\n");
	#endif
	
	
	printf( "nDeme..... %d\n", p.nDeme);
	printf( "nTry...... %d\n", p.nTry);
	printf( "nGen...... %d\n", p.nGen);
	
	printf( "nInd1_2_T. ");
	for (i =0 ; i < p.nDeme; i++)
	{
		printf("%d ", p.popsize[i]);
	}
	printf("%d \n", p.nIndT);
	
	
	printf( "N_S_Total. %d %d %d\n", p.nNon, p.nSyn, p.nSeq);
	
	printf( "mig....... %d\n", p.mig);
	printf( "period.... %d\n", p.period);

	
	printf( "Nu........ %lf\n", p.Nu);
	printf( "recomRate. %lf \n", p.recomRate);

	printf( "sMax...... %lf\n", p.sMax);
	printf( "epsilon... %lf\n", p.epsilon);
	printf( "C......... %lf\n", p.C);
	printf( "sigma_s... %lf\n", p.sigma_s);
	printf( "del_range. %lf\n", p.del_range);
	
	printf( "seed...... %ld\n", p.seed);
	
	printf("==================================\n\n");
	
	
	
	//print parameters to outFile

	#if defined(_SBC_)
	fprintf(g_outfile, " << SEEDBANK - CONTROL    MODE >> \n");
	#elif defined(_SBE_)
	fprintf(g_outfile, " << SEEDBANK - EXPERIMNET MODE >> \n");
	#elif defined(_RSC_)
	fprintf(g_outfile, " << REFUGE   - CONTROL    MODE >> \n");
	#elif defined(_RSE_)
	fprintf(g_outfile, " << REFUGE   - EXPERIMNET MODE >> \n");
	#else 
	fprintf(g_outfile, "++++++++++++++++++++++++++++++\n+[Warning] Please compile with option+\n++++++++++++++++++++++++++++++\n");
	#endif
	
	
	fprintf(g_outfile, "nDeme %d\n", p.nDeme);
	fprintf(g_outfile, "nTry %d\n", p.nTry);
	fprintf(g_outfile, "nGen %d\n", p.nGen);
	
	//fprintf(g_outfile, "nInd1_2_Total %d %d %d\n", p.nInd1,  p.nInd1, p.nIndT);
	fprintf(g_outfile, "nInd1_2_T. ");
	for (i =0 ; i < p.nDeme; i++)
	{
		fprintf(g_outfile,"%d ", p.popsize[i]);
	}
	fprintf(g_outfile,"%d \n", p.nIndT);

	fprintf(g_outfile, "N_S_Total %d %d %d\n", p.nNon, p.nSyn, p.nSeq);
	
	fprintf(g_outfile, "mig %d\n", p.mig);
	fprintf(g_outfile, "period %d\n", p.period);

	
	fprintf(g_outfile, "Nu %lf\n", p.Nu);
	fprintf(g_outfile, "recomRate %lf\n", p.recomRate);

	
	fprintf(g_outfile, "sMax %lf\n", p.sMax);
	fprintf(g_outfile, "epsilon %lf\n", p.epsilon);
	fprintf(g_outfile, "C %lf\n", p.C);
	fprintf(g_outfile, "sigma_s %lf\n", p.sigma_s);
	fprintf(g_outfile, "del_range. %lf\n", p.del_range);
	
	fprintf(g_outfile, "seed %ld\n", p.seed);
	
	fprintf(g_outfile,"==================================\n\n");


	return p;
}


void AssignNonsynSite(unsigned int* template_nonsys, para_t p)
{
	//make n sites to be nonsynonymous sites
	//mark nonsynonymous site as 1
	//return template
	
	int i, j, count = 0, bit;
	int nNonsyn, rSite;
	int r_site, r_block;
	
	unsigned int seg;

	int nSeq = p.nSeq;
	
	memset(template_nonsys, 0, sizeof(unsigned int) * p.nSeq_block);
	
	nNonsyn = p.nNon;
	while (count < nNonsyn)
	{
		rSite = (int) (ran1(&gseed)*nSeq);	
		g_rSite = rSite;
		
		r_block = rSite/_block_;
		r_site = rSite%_block_;
		
		seg = template_nonsys[r_block];
		bit = (seg &  1 << (r_site))?1:0;
		
		
		
																		#if PRINT_TML 
																		printf("\nnonNum = %d block %2d site %2d  (%d) >>%d<<\n" , count, r_block, r_site,rSite, bit);
																		for (j = 0 ; j < p.nSeq_block; j++)
																		{
																			PrintBinary_asterik(template_nonsys[j]);
																			printf(" ");
																		}printf("\n");
																		
																		if (r_block == 1)
																		{
																			printf("r_block = 1\n");
																			PrintBinary_asterik(0);
																			printf(" ");
																			PrintBinary_asterik(1 << (r_site));
																			printf("\n");
																		}
																		else if (r_block == 0)
																		{
																			printf("r_block = 0\n");
																			PrintBinary_asterik(1 << (r_site));
																			printf(" ");
																			PrintBinary_asterik(0);
																			printf("\n");
																			
																			
																		}
																		#endif
		
		if ( bit == 0 )
		{
			seg += ( 1 << (r_site) );
			count++;
		}

		template_nonsys[r_block] = seg;

		
																		#if PRINT_TML 
														
																		for (j = 0 ; j < p.nSeq_block; j++)
																		{
																			PrintBinary_asterik(template_nonsys[j]);
																			printf(" ");
																		}
																		#endif
	
	}
	
}
	
void AssignDeltaEachNon(double* delta, para_t p)
{
	int i;
	
	memset(delta, 0, sizeof(double)* p.nNon);
	for (i = 0; i<p.nNon; i++)
	{
		#if DELTA_FIXED
		delta[i] = _DELTA_ * pow(-1, i);
		#else
		delta[i] = p.del_range * gasdev(&gseed);
		#endif
	}
												#if PRINT_TML
												printf("delta : " );
												for (i = 0; i<p.nNon; i++)
												{
													printf("%lf ", delta[i]);
												}
												printf("\n");
												#endif	
}

void AssignDeltaEachNon_Exp(double* delta, para_t p)
{
	int i;
	
	memset(delta, 0, sizeof(double)* p.nNon);
	for (i = 0; i<p.nNon; i++)
	{
		#if DELTA_FIXED
		delta[i] = _DELTA_ * pow(-1, i);
		#else
		delta[i] = p.del_range * expdev(&gseed);
		#endif
	}
												#if PRINT_TML
												printf("delta : " );
												for (i = 0; i<p.nNon; i++)
												{
													printf("%lf ", delta[i]);
												}
												printf("\n");
												#endif	
}




void AssignSelOptima(double** selOpt, int mode, para_t p)
{
	int i, j;
	int r = gseed%p.period;
	double s_t;
	int half_period = p.period/2;
	
	
	for (i =0;i<p.nDeme;i++)
	{
		
			memset(selOpt[i], 0, sizeof(double) * p.period);
		
	}
	
	
	//printf(" p.sMax = %lf\n",  p.sMax);
	if (mode == 1 )	// binary
	{
		
		for (i = 0; i < p.nDeme; i++)
		{
			for (j = 0; j < p.period; j++)
			{
				if (i ==0)
				{
					
					if ( j < half_period )
						selOpt[i][j] =  p.sMax;
					else
						selOpt[i][j] = - p.sMax;
				}
				else
				{
					if ( j < half_period )
						selOpt[i][j] = p.C * p.sMax;
					else
						selOpt[i][j] = - (p.C * p.sMax);
				}
		
			}
		}
	}
	else if (mode == 2) // sine fuction
	{
		for (i = 0; i < p.nDeme; i++)
		{
			for (j = 0; j < p.period; j++)
			{
				s_t =  p.sMax * sin(2*PI*((double)(i+r)/p.period));	
				if (i ==0)
				{
					selOpt[i][j] =  s_t;
				}
				else
				{
					selOpt[i][j] =  p.C  * s_t;
				}
				
			}		
				
		}
	}
											#if PRINT_TML
												//for (i = 0; i < p.nDeme; i++)
												//{
												//	printf("subPop%d ------------ \n", i);

													for (j = 0; j < p.period; j++)
													{
														printf("%dth  : %lf\t%lf\n", j, selOpt[0][j], selOpt[1][j]);
													}
													//printf("---------------------\n");
												//}
											#endif
}

void DisturbSelOpt(double** selOpt,double** selOpt_origianl, double epsilon, int period, int nDeme)
{
	//only works for Gul model -> add same Epsilon value to sub1 and sub2
	int i,j;
	int half_period = period/2;
	double ranE1, ranE2;
	
	ranE1 = (ran1(&gseed) * 2 * epsilon) - epsilon;	//random from (-epsilon, epsilon
	ranE2 = (ran1(&gseed) * 2 * epsilon) - epsilon;

	
	for (i = 0; i < nDeme ; i++)
	{
		for (j = 0; j < period; j++)
		{
			
			if ( j < half_period )
				selOpt[i][j] = selOpt_origianl[i][j] + ranE1;
			else
				selOpt[i][j] = selOpt_origianl[i][j] + ranE2;
			
	
		}
	}
	
	
	
	
											#if PRINT_TML
												//for (i = 0; i < p.nDeme; i++)
												//{
													printf("=====during this period, ===== \n" );
													printf(" ranE1 : %lf, ranE2 : %lf  \n", ranE1, ranE2);

													for (j = 0; j < period; j++)
													{
														printf("%dth  : %lf\t%lf\n", j, selOpt[0][j], selOpt[1][j]);
													}
													printf("---------------------\n\n");
												//}
											#endif
}


struct Population* allocPop(para_t p)
{
	#if PRINT_WHERE
	printf("\n\n==========allocPop==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========allocPop==========\n");
	#endif
	
	int dm, i;
	
	
	struct Population* pop = (struct Population*) malloc ( sizeof(struct Population) * p.nDeme);
	
	for (dm=0; dm<p.nDeme; dm++)
	{
		pop[dm].indiv = (struct Individual*)malloc(sizeof(struct Individual) * p.popsize[dm]);
		
		for(i = 0; i<p.popsize[dm]; i++)
		{
			pop[dm].indiv[i].seq = (unsigned int*) malloc(p.nSeq_block*sizeof(unsigned int));
			memset(pop[dm].indiv[i].seq,  0, sizeof(unsigned int)*p.nSeq_block);
		}
	}
	return pop;
}



double FitOfIndiv(struct Individual* ind, unsigned int* nonsyn, double* delta, double pheno_wild, double selOpt, para_t p)
{
	double pheno, fit;
	double sigma_s = p.sigma_s;
	int* seq = ind->seq;
		
	pheno = getPhenoFromSeq(seq, nonsyn, delta, pheno_wild, p);
	//fit = GetFitnessFromPheno(pheno, selOpt, gen);
	
	
	//formola from Svardal2011
	//Gaussian stabilizing selection
	fit = exp(-(pheno - selOpt)*(pheno - selOpt) / (2.0 * sigma_s * sigma_s));
	#if PRINT_FIT
	printf("z = %lf / selOpt = %lf -> exp(%lf) = %lf ", pheno, selOpt, (-(pheno - selOpt)*(pheno - selOpt) / (2.0 * sigma_s * sigma_s)), fit);
	#endif
	
	return fit;
}


double FitOfIndiv_noStorage(struct Individual* ind, unsigned int* nonsyn, double* delta, double pheno_wild, double selOpt, para_t p)
{
	double pheno, fit;
	double sigma_s = p.sigma_s;
	int* seq = ind->seq;
		
	pheno = getPhenoFromSeq(seq, nonsyn, delta, pheno_wild, p);
	//fit = GetFitnessFromPheno(pheno, selOpt, gen);
	
	
	selOpt = 0.5 * selOpt;
	//formola from Svardal2011
	//Gaussian stabilizing selection
	fit =  exp( - (((pheno - selOpt) * (pheno - selOpt)) / (2.0 * sigma_s * sigma_s)));
	#if PRINT_FIT
	printf("z = %lf / selOpt = %lf -> exp(%lf) = %lf ", pheno, selOpt, (-(pheno - selOpt)*(pheno - selOpt) / (2.0 * sigma_s * sigma_s)), fit);
	#endif
	
	return fit;
}


/*
void Recombination(struct Individual* ind1, struct Individual* ind2, int recombinantNum, struct Individual* recombinant, para_t p)
{
	//In one segment(=seq_block), recombination rate is constant
	//recombination can occur anywhere, even inside the segment
	#if PRINT_RECOMB
	printf("\n\n==========Recombination==========\n\n");
	#endif
	
	int n, s;
	unsigned int seg1, seg2, tmp_seg1, tmp_seg2, template_recom;
	int nRecom, rSite;
	
	//
	double recomRate = p.recomRate;
	

	for (s = 0 ; s < p.nSeq_block; s++)
	{
		nRecom = bnldev(recomRate, _block_, &gseed);	//number of recombination in a block(=segment)
		seg1 = ind1->seq[s];
		seg2 = ind2->seq[s];
											#if PRINT_RECOMB
											printf("in block %d, nRecom = %d	\n", s, nRecom);		
											printf(" before Recombination			after Recombination\n");
											#endif
		for (n = 0; n < nRecom ; n++)
		{

			template_recom = 0;
			rSite = (int) (ran1(&gseed)*_block_);	
			template_recom = ((1 << rSite) - 1);
		
			tmp_seg1 = ((~template_recom) & seg1) + ((template_recom) & seg2);
			tmp_seg2 = ((~template_recom) & seg2) + ((template_recom) & seg1);
			
											#if PRINT_RECOMB
											PrintBinary(seg1);
											printf(" -%2d- > ", rSite);
											PrintBinary(tmp_seg1);
											printf("\n");
											
											template_recom = (1 << (rSite)) - 1;
											PrintBinary_asterik(template_recom);
											printf(" -  - > ");
											PrintBinary_asterik(template_recom);
											printf("\n");
											
											PrintBinary(seg2);
											printf(" -  - > ");
											PrintBinary(tmp_seg2);
											printf("\n");
											printf("\n");
											#endif

			
			seg1 = tmp_seg1;
			seg2 = tmp_seg2;
			
		}
	}
				
}

*/


//void Recombination(struct Individual* ind1, struct Individual* ind2, para_t p)
void Recombination(struct Individual* ind1, struct Individual* ind2, int recombinantNum, struct Individual* recombinant, para_t p)
{
	//In one segment(=seq_block), recombination rate is constant
	//recombination can occur anywhere, even inside the segment
	
	#if PRINT_RECOMB
	printf("\n\n==========Recombination==========\n\n");
	#endif
	
	unsigned int* template_recom;
	unsigned int template, tmp_seg1, tmp_seg2, seg1, seg2, segR;
	
	double recomRate = p.recomRate;
	
	int rSite, nRecom, arr=0;

	int i, n, s;
	
	if (recomRate == 0)			//for no recom simulation
	{
		if (recombinantNum == 0)
		{
			for (s = 0; s < p.nSeq_block; s++)
			{
				recombinant->seq[s] = ind1->seq[s];
			}
		}
		else
		{
			for (s = 0; s < p.nSeq_block; s++)
			{
				recombinant->seq[s] = ind2->seq[s];
			}
		}
		return;
	}
	
	  
	nRecom = poidev(recomRate * p.nSeq , &gseed);			
													#if PRINT_RECOMB
													printf("nRecom = %d	\n", nRecom);
													printf(" before after \n");
													#endif
	
	// tem array 할당하기

	template_recom = (unsigned int*) malloc (p.nSeq_block * sizeof(unsigned int));
	for (i = 0; i < p.nSeq_block; i++)
	{
		template_recom[i] = 0;
	}
	
	#if PRINT_RECOMB
	printf("Array malloc\n");
	#endif
  
	
	//making a template for recombination
	for (n = 0; n < nRecom; n++)
	{
		rSite = (int) (ran1(&gseed) * p.nSeq);

													#if PRINT_RECOMB
													printf("| rSite %d	block %d	site %d |\n", rSite, rSite/_block_, rSite%_block_);		
													
													for (s = 0; s < p.nSeq_block; s++)
													{
														PrintBinary(template_recom[s]);
														printf(" ");
													}
													printf("\n");
													
																				
													
													#endif
													
													
		for (s = 0; s < p.nSeq_block; s++)
		{
			if ( s > (rSite/_block_) )										// blocks after the block that recombination happens
			{
				template = ( 1 << (_block_)) - 1 ;							//template 		= 11111
				template_recom[s] = ( template_recom[s] ^ template );		//template[s]	= 00000 -> 11111
																			// or 		  	= 01010 -> 01010
																			
													#if PRINT_RECOMB
													PrintBinary_asterik(template);
													printf(" ");
													#endif
			}
				
			else if ( s == (rSite/_block_) )								// the block that recombination happens
			{
				template = ( 1 << ((rSite%_block_)+1) ) -1;					//template 		= 00111 if rSite%block = 3
																			// or 			= 11111 if rSite%block = 5
				template_recom[s] = ( template_recom[s] ^ template );		//template[s]	= 00000 -> 00111
																			// or 		  	= 01010 -> 01101
													#if PRINT_RECOMB
													PrintBinary_asterik(template);
													printf(" ");
													#endif
			}
			
													#if PRINT_RECOMB
													else
													{	
														PrintBinary_asterik(0);
														printf(" ");
													}
													#endif
													
																			// blocks before the block that recombination happens			
		}

		
		
													#if PRINT_RECOMB
													printf("\n");
													for (s = 0; s < p.nSeq_block; s++)
													{
														PrintBinary(template_recom[s]);
														printf(" ");
													}
													printf("\n............................................................................\n");
													#endif
	}
	
											
	  												
	//make recombinants using template
	//seg1 = ind1->seq[s];				//여기 에러였던듯!!! s에 값이 없어서 할당이 안되었었음!!
	//seg2 = ind2->seq[s];
													#if PRINT_RECOMB
													for (s = 0; s < p.nSeq_block; s++)
													{
														PrintBinary(ind1->seq[s]);
														printf(" ");
													}
													printf(" original 1\n");
													for (s = 0; s < p.nSeq_block; s++)
													{
														PrintBinary_asterik(template_recom[s]);
														printf(" ");
													}
													printf("   template\n");
													for (s = 0; s < p.nSeq_block; s++)
													{
														PrintBinary(ind2->seq[s]);
														printf(" ");
													}
													printf(" original 2\n~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
													#endif
	for (s = 0; s < p.nSeq_block; s++) 
	{
		//seg1 = ind1->seq[s];
		//seg2 = ind2->seq[s];
		//segR = recombinant->seq[s];
		
		tmp_seg1 = ((~template_recom[s]) & ind1->seq[s]) + ((template_recom[s]) & ind2->seq[s]);
		tmp_seg2 = ((~template_recom[s]) & ind2->seq[s]) + ((template_recom[s]) & ind1->seq[s]);
	
		if (recombinantNum == 0)
			recombinant->seq[s] = tmp_seg1;
		else	
			recombinant->seq[s] = tmp_seg2;
		//seg2 = tmp_seg2;
		
	}
													#if PRINT_RECOMB
													
													for (s = 0; s < p.nSeq_block; s++)
													{
														PrintBinary(((~template_recom[s]) & ind1->seq[s]) + ((template_recom[s]) & ind2->seq[s]));
														printf(" ");
													}
													if (recombinantNum == 0) 
														printf(" recombin 1     V\n");
													else
														printf(" recombin 1\n");
													for (s = 0; s < p.nSeq_block; s++)
													{
														PrintBinary(((~template_recom[s]) & ind2->seq[s]) + ((template_recom[s]) & ind1->seq[s]));
														printf(" ");
													}
													if (recombinantNum == 1)
														printf(" recombin 2     V\n---------------------------\n");
													else
														printf(" recombin 2\n---------------------------\n");
													
													for (s = 0; s < p.nSeq_block; s++)
													{
														PrintBinary(recombinant->seq[s]);
														printf(" ");
													}printf(" recombinant**       \n");
													#endif
	
	
	free(template_recom); 
	#if PRINT_RECOMB
			
	
		
		printf("===========Array free===========\n");
	
	#endif	
	  
	//인쇄해서 잘 돌아가는지 확인해야함 + method에 정리해야 함
}

/*
void Recombinant(struct Individual* ind1, struct Individual* ind2, struct Individual* recombinant, para_t p)
{
	int i;
	for ( i = 0; i < popSize; i++)
	{
		recombinant[i] = i;
	}
	
	
	return;
	
}
	*/


	

struct Individual CopyIndividual(struct Individual dest, struct Individual sourc, para_t p)
{
	memcpy(dest.seq, sourc.seq, (sizeof(int))*p.nSeq_block);
	//indiv0 : destination / indiv1 : source
	return dest;
	
}

void Migration(struct Population* pop, para_t p)
{
	#if PRINT_WHERE
	printf("\n\n==========Migration==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========Migration==========\n");
	#endif
	
	struct Individual tmp;
	tmp.seq = (unsigned int*) malloc(sizeof(unsigned int)*p.nSeq_block);
	
	int i;
	
	if (p.nDeme == 2)
	{
		for (i=0; i < p.mig; i++)
		{
			tmp = CopyIndividual(tmp, pop[0].indiv[i], p);
			pop[0].indiv[i] = CopyIndividual(pop[0].indiv[i], pop[1].indiv[i], p);
			pop[1].indiv[i] = CopyIndividual(pop[1].indiv[i] ,tmp ,p);
		}
	}
	free(tmp.seq);
}





/*
struct Population* SelectionAndReproducation(struct Population* currentPop, struct Population* pop1, struct Population* pop2,unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen, para_t p)
{
	#if PRINT_WHERE
	printf("\n\n==========SelectionAndReproducation==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========SelectionAndReproducation==========\n");
	#endif
	
	
	struct Population* parentPop;
	struct Population* offspringPop;
	int rInd1, rInd2;
	int i,j, count;
	int gen_period = gen % (p.period);
	double fit1, fit2;
	double oldpheno1=-100000,oldpheno2=-100000,  nowpheno;
	int phenoCount1=0, phenoCount2=0;
	int* seq;
	
	if (currentPop == pop1)								//pop1 = parentPop, pop2 = offspringPop
	{
		parentPop = pop1;
		offspringPop = pop2;
	}
	else if (currentPop == pop2)						//pop2 = parentPop, pop1 = offspringPop	
	{
		parentPop = pop2;
		offspringPop = pop1;		
	}
	
	
	
	for (i =0 ; i < p.nDeme; i++)
	{
		count = 0;
		while (count < p.nInd1)
		{
																						
			rInd1 = (int) (p.nInd1 * ran1(&gseed));
			rInd2 = (int) (p.nInd1 * ran1(&gseed));
			//Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), p);
			
			
																						#if PRINT_SELECTION
																						printf("subpop %d -----------------------------------------------------%d\n",i,  count);
																						
																						printf("Ind1(parent : %5d) :", rInd1);
																						for (j=0;j<p.nSeq_block;j++){
																							PrintBinary(parentPop[i].indiv[rInd1].seq[j]);
																							printf(" ");}
																						printf("\nInd2(parent : %5d) :", rInd2);
																						for (j=0;j<p.nSeq_block;j++){
																							PrintBinary(parentPop[i].indiv[rInd2].seq[j]);
																							printf(" ");}
																						printf("\n");
																					//	printf("> %lf		%lf \n	", PhenoOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, selOpt[i][gen_period], 1,p), PhenoOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, selOpt[i][gen_period], 1,p));
																						
																						#endif
		
			fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
			fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
			
			
			
			if(ran1(&gseed) < fit1 )
			//	ran1 < relative fitness
			// = ran1 * mostFit < indFit
			// = ran1 * 1 < indFit  (since mostFit = 1, see GetFitnessFromPheno)
			{
				#if TRJ
					seq = parentPop[i].indiv[rInd1].seq;
					nowpheno = getPhenoFromSeq(seq, nonsyn, delta, pheno_wild,p);
					if (oldpheno1 == nowpheno)
					{
						phenoCount1++;
						
					}
					else if(oldpheno2 == nowpheno)
					{
						phenoCount2++;
					}
					
					else
					{	
						if (oldpheno2!= -100000)
							fprintf(g_phenofile, "%lf_%d ", oldpheno2 , phenoCount2);
						
						oldpheno2 = oldpheno1;
						oldpheno1 = nowpheno;
						phenoCount2 = phenoCount1;
						phenoCount1 = 1;
					}
					
					
				#endif
				
				offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[rInd1], p);
				count +=1;
																						#if PRINT_SELECTION
																						printf("	O" );
																						
																						if (count == p.nInd1){
																							printf("\noffspringPop\n");
																							PrintPop(offspringPop, p);
																							break;}
																							
																						#else
				if (count == p.nInd1) break;
																						#endif	
			}
			
			if(ran1(&gseed) < fit2 )
			{
				#if TRJ
					seq = parentPop[i].indiv[rInd2].seq;
					nowpheno = getPhenoFromSeq(seq, nonsyn, delta, pheno_wild,p);
					if (oldpheno1 == nowpheno)
					{
						phenoCount1++;
						
					}
					else if(oldpheno2 == nowpheno)
					{
						phenoCount2++;
					}
					
					else
					{	
						if (oldpheno2!= -100000)
							fprintf(g_phenofile, "%lf_%d ", oldpheno2 , phenoCount2);
						
						oldpheno2 = oldpheno1;
						oldpheno1 = nowpheno;
						phenoCount2 = phenoCount1;
						phenoCount1 = 1;
					}
					
					
				#endif
				offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[rInd2], p);
				count +=1;
																						#if PRINT_SELECTION
																						printf("	O" );
																					
																						if (count == p.nInd1){
																							printf("\noffspringPop\n");
																							PrintPop(offspringPop, p);
																							break;}
																						#else
				if (count == p.nInd1) break;
																						#endif	
			}
			
			
																						#if PRINT_SELECTION
																						printf("\noffspringPop\n");
																						PrintPop(offspringPop, p);
																						#endif
			
			
			
		}
				
			
	}

	#if TRJ
			if (oldpheno2!= -100000)
				fprintf(g_phenofile, "%lf_%d ", oldpheno2 , phenoCount2);
			if (oldpheno1!= -100000)
				fprintf(g_phenofile, "%lf_%d ", oldpheno1 , phenoCount1);
			#endif	
	return offspringPop;
}

struct Population* SelectionAndReproducation_Guli(struct Population* currentPop, struct Population* pop1, struct Population* pop2,unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen, para_t p)
{
	#if PRINT_WHERE
	printf("\n\n==========SelectionAndReproducation==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========SelectionAndReproducation==========\n");
	#endif
	
	
	struct Population* parentPop;
	struct Population* offspringPop;
	int rInd1, rInd2;
	int i,j, count;
	int gen_period = gen % (p.period);
	double fit1, fit2;
	double oldpheno1=-100000,oldpheno2=-100000,  nowpheno;
	int phenoCount1=0, phenoCount2=0;
	int* seq;
	
	if (currentPop == pop1)								//pop1 = parentPop, pop2 = offspringPop
	{
		parentPop = pop1;
		offspringPop = pop2;
	}
	else if (currentPop == pop2)						//pop2 = parentPop, pop1 = offspringPop	
	{
		parentPop = pop2;
		offspringPop = pop1;		
	}
	                                     
	
	
	for (i =0 ; i < p.nDeme; i++)
	{
		count = 0;
		while (count < p.nInd1)
		{
																						
			rInd1 = (int) (p.nInd1 * ran1(&gseed));
			rInd2 = (int) (p.nInd1 * ran1(&gseed));
			//Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), p);
			
			
																						#if PRINT_SELECTION
																						printf("subpop %d -----------------------------------------------------%d\n",i,  count);
																						
																						printf("Ind1(parent : %5d) :", rInd1);
																						for (j=0;j<p.nSeq_block;j++){
																							PrintBinary(parentPop[i].indiv[rInd1].seq[j]);
																							printf(" ");}
																						printf("\nInd2(parent : %5d) :", rInd2);
																						for (j=0;j<p.nSeq_block;j++){
																							PrintBinary(parentPop[i].indiv[rInd2].seq[j]);
																							printf(" ");}
																						printf("\n");
																					//	printf("> %lf		%lf \n	", PhenoOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, selOpt[i][gen_period], 1,p), PhenoOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, selOpt[i][gen_period], 1,p));
																						
																						#endif
		
			if (i == 0)
			{
				//fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				//fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				
				#if noStorageMode
				//fit1 = FitOfIndiv_noStorage(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				//fit2 = FitOfIndiv_noStorage(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				
				#else
				fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				#endif
			}
			else if (i == 1)
			{
				#if noStorageMode
				//fit1 = FitOfIndiv_noStorage(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				//fit2 = FitOfIndiv_noStorage(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				#else
				fit1 = 1;
				fit2 = 1;
				#endif
			}
			
			
			if(ran1(&gseed) < fit1 )
			//	ran1 < relative fitness
			// = ran1 * mostFit < indFit
			// = ran1 * 1 < indFit  (since mostFit = 1, see GetFitnessFromPheno)
			{
				#if TRJ
					seq = parentPop[i].indiv[rInd1].seq;
					nowpheno = getPhenoFromSeq(seq, nonsyn, delta, pheno_wild,p);
					if (oldpheno1 == nowpheno)
					{
						phenoCount1++;
						
					}
					else if(oldpheno2 == nowpheno)
					{
						phenoCount2++;
					}
					
					else
					{	
						if (oldpheno2!= -100000)
							fprintf(g_phenofile, "%lf_%d ", oldpheno2 , phenoCount2);
						
						oldpheno2 = oldpheno1;
						oldpheno1 = nowpheno;
						phenoCount2 = phenoCount1;
						phenoCount1 = 1;
					}
					
					
				#endif
				
				offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[rInd1], p);
				count +=1;
																						#if PRINT_SELECTION
																						printf("	O" );
																						
																						if (count == p.nInd1){
																							printf("\noffspringPop\n");
																							PrintPop(offspringPop, p);
																							break;}
																							
																						#else
				if (count == p.nInd1) break;
																						#endif	
			}
			
			if(ran1(&gseed) < fit2 )
			{
				#if TRJ
					seq = parentPop[i].indiv[rInd2].seq;
					nowpheno = getPhenoFromSeq(seq, nonsyn, delta, pheno_wild,p);
					if (oldpheno1 == nowpheno)
					{
						phenoCount1++;
						
					}
					else if(oldpheno2 == nowpheno)
					{
						phenoCount2++;
					}
					
					else
					{	
						if (oldpheno2!= -100000)
							fprintf(g_phenofile, "%lf_%d ", oldpheno2 , phenoCount2);
						
						oldpheno2 = oldpheno1;
						oldpheno1 = nowpheno;
						phenoCount2 = phenoCount1;
						phenoCount1 = 1;
					}
					
					
				#endif
				offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[rInd2], p);
				count +=1;
																						#if PRINT_SELECTION
																						printf("	O" );
																					
																						if (count == p.nInd1){
																							printf("\noffspringPop\n");
																							PrintPop(offspringPop, p);
																							break;}
																						#else
				if (count == p.nInd1) break;
																						#endif	
			}
			
			
																						#if PRINT_SELECTION
																						printf("\noffspringPop\n");
																						PrintPop(offspringPop, p);
																						#endif
			
			
			
		}
				
			
	}

	#if TRJ
			if (oldpheno2!= -100000)
				fprintf(g_phenofile, "%lf_%d ", oldpheno2 , phenoCount2);
			if (oldpheno1!= -100000)
				fprintf(g_phenofile, "%lf_%d ", oldpheno1 , phenoCount1);
			#endif	
	return offspringPop;
}




struct Population* SelectionAndReproducation_Guli_SeedBank(struct Population* currentPop, struct Population* pop1, struct Population* pop2,unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen, para_t p)
{
	#if PRINT_WHERE
	printf("\n\n==========SelectionAndReproducation==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========SelectionAndReproducation==========\n");
	#endif
	
	
	struct Population* parentPop;
	struct Population* offspringPop;
	int rInd1, rInd2, randInd;
	int i,j, count, k, tmp;
	int gen_period = gen % (p.period);
	double fit1, fit2;
	double oldpheno1=-100000,oldpheno2=-100000,  nowpheno;
	int phenoCount1=0, phenoCount2=0;
	int* seq;
	int* ranArray;
	
	if (currentPop == pop1)								//pop1 = parentPop, pop2 = offspringPop
	{
		parentPop = pop1;
		offspringPop = pop2;
	}
	else if (currentPop == pop2)						//pop2 = parentPop, pop1 = offspringPop	
	{
		parentPop = pop2;
		offspringPop = pop1;		
	}
	
	ranArray = (int*) malloc (p.nInd1 * sizeof(int) );		//for random sorting of second subpopulation 
	for ( i = 0; i < p.nInd1; i++)
	{
		ranArray[i] = i;
	}
	
	
	//First Subpop---------------------------------
	
	
	for (i =0 ; i < p.nDeme; i++)
	{
		count = 0;
		while (count < p.nInd1)
		{
						
			if (i == 0)
			{
				rInd1 = (int) (p.nInd1 * ran1(&gseed));
				rInd2 = (int) (p.nInd1 * ran1(&gseed));
				//Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), p);
				
				
																							#if PRINT_SELECTION
																							printf("\nsubpop %d -----------------------------------------------------%d\n",i,  count);
																							
																							printf("Ind1(parent : %5d) :", rInd1);
																							for (j=0;j<p.nSeq_block;j++){
																								PrintBinary(parentPop[i].indiv[rInd1].seq[j]);
																								printf(" ");}
																							printf("\nInd2(parent : %5d) :", rInd2);
																							for (j=0;j<p.nSeq_block;j++){
																								PrintBinary(parentPop[i].indiv[rInd2].seq[j]);
																								printf(" ");}
																							printf("\n");
																						//	printf("> %lf		%lf \n	", PhenoOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, selOpt[i][gen_period], 1,p), PhenoOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, selOpt[i][gen_period], 1,p));
																							
																							#endif
			
				
				

					
				#if noStorageMode
				fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
					
				#else
				fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
				#endif
					
					
					
					
				if(ran1(&gseed) < fit1 )
				//	ran1 < relative fitness
				// = ran1 * mostFit < indFit
				// = ran1 * 1 < indFit  (since mostFit = 1, see GetFitnessFromPheno)
				{
					#if TRJ
						seq = parentPop[i].indiv[rInd1].seq;
						nowpheno = getPhenoFromSeq(seq, nonsyn, delta, pheno_wild,p);
						if (oldpheno1 == nowpheno)
						{
							phenoCount1++;
							
						}
						else if(oldpheno2 == nowpheno)
						{
							phenoCount2++;
						}
						
						else
						{	
							if (oldpheno2!= -100000)
								fprintf(g_phenofile, "%lf_%d ", oldpheno2 , phenoCount2);
							
							oldpheno2 = oldpheno1;
							oldpheno1 = nowpheno;
							phenoCount2 = phenoCount1;
							phenoCount1 = 1;
						}
						
						
					#endif
					
					offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[rInd1], p);
					count +=1;
																							#if PRINT_SELECTION
																							printf("	O" );
																							
																							if (count == p.nInd1){
																								printf("\noffspringPop\n");
																								PrintPop(offspringPop, p);
																								break;}
																								
																							#else
					if (count == p.nInd1) break;
																							#endif	
				}
				
				if(ran1(&gseed) < fit2 )
				{
					#if TRJ
						seq = parentPop[i].indiv[rInd2].seq;
						nowpheno = getPhenoFromSeq(seq, nonsyn, delta, pheno_wild,p);
						if (oldpheno1 == nowpheno)
						{
							phenoCount1++;
							
						}
						else if(oldpheno2 == nowpheno)
						{
							phenoCount2++;
						}
						
						else
						{	
							if (oldpheno2!= -100000)
								fprintf(g_phenofile, "%lf_%d ", oldpheno2 , phenoCount2);
							
							oldpheno2 = oldpheno1;
							oldpheno1 = nowpheno;
							phenoCount2 = phenoCount1;
							phenoCount1 = 1;
						}
						
						
					#endif
					offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[rInd2], p);
					count +=1;
																							#if PRINT_SELECTION
																							printf("	O" );
																						
																							if (count == p.nInd1){
																								printf("\noffspringPop\n");
																								PrintPop(offspringPop, p);
																								break;}
																							#else
					if (count == p.nInd1) break;
																							#endif	
				}
				
			
			

			}
			else if (i == 1)
			{
																						#if PRINT_SELECTION
																						printf("subpop %d -----------------------------------------------------%d\n",i,  count);
																						#endif
				for (k = 0; k < p.nInd1; k++)
				{
					randInd = ran1(&gseed) * (p.nInd1 - k) + k ; 
					
					tmp = ranArray[k];
					ranArray[k] = ranArray[randInd];
					ranArray[randInd] = tmp;
																						#if PRINT_SELECTION
																						printf("%d ", ranArray[k]);
																						#endif
					#if TRJ
					seq = parentPop[i].indiv[ranArray[k]].seq;
					nowpheno = getPhenoFromSeq(seq, nonsyn, delta, pheno_wild,p);
					if (oldpheno1 == nowpheno)
					{
						phenoCount1++;
						
					}
					else if(oldpheno2 == nowpheno)
					{
						phenoCount2++;
					}
					
					else
					{	
						if (oldpheno2!= -100000)
							fprintf(g_phenofile, "%lf_%d ", oldpheno2 , phenoCount2);
						
						oldpheno2 = oldpheno1;
						oldpheno1 = nowpheno;
						phenoCount2 = phenoCount1;
						phenoCount1 = 1;
					}
					
						
					#endif															
																						
					
					
					
					offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[ranArray[k]], p);
					count +=1;
				}
																						
																						#if PRINT_SELECTION
																						printf("\n");
																				
																						#endif
				if (count == p.nInd1) break;
					
				
			}
			
			
			
			
	
			
			
		}
				
			
	}

			#if TRJ
			if (oldpheno2!= -100000)
				fprintf(g_phenofile, "%lf_%d ", oldpheno2 , phenoCount2);
			if (oldpheno1!= -100000)
				fprintf(g_phenofile, "%lf_%d ", oldpheno1 , phenoCount1);
			#endif	
			
	free(ranArray);
	
	return offspringPop;
}


*/



struct Population* SelectionAndReproducation_Guli_SeedBank(struct Population* currentPop, struct Population* pop1, struct Population* pop2,unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen, list* haploTable, para_t p)
{
	/* 	Based on Wright-Fisher model
		
		In subpopulation 1, 
		  randomly select two individuals and get fitness
		  if ranNum > fitness, then the individual undergoes recombination with the other indiv - for both indiv
		  offspring will be recombinant of the individual and haploType is counted
		  repeat until count == popSize ( the population size of subPop1)
		
		In subpopulation 2 (SeedBank), 
		  the order of individual is randomly reassigned in offspring population using for loop
		  while copying indivs to offpop, haploType is counted
	*/	  
	
	
	#if PRINT_WHERE
	printf("\n\n==========SelectionAndReproducation==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========SelectionAndReproducation==========\n");
	#endif
	
	
	struct Population* parentPop;
	struct Population* offspringPop;
	int rInd1, rInd2, randInd;
	int i,j, count, k, tmp;
	int gen_period = gen % (p.period);
	double fit1, fit2;
	double oldpheno1=-100000,oldpheno2=-100000,  nowpheno;
	int phenoCount1=0, phenoCount2=0;
	int* seq;
	int* ranArray;
	int popSize;
	
	
	if (currentPop == pop1)								//pop1 = parentPop, pop2 = offspringPop
	{
		parentPop = pop1;
		offspringPop = pop2;
	}
	else if (currentPop == pop2)						//pop2 = parentPop, pop1 = offspringPop	
	{
		parentPop = pop2;
		offspringPop = pop1;		
	}
	
	
	
	struct Individual recombinant;
	recombinant.seq = (unsigned int*) malloc (p.nSeq_block * sizeof(unsigned int));
	
	
	
	//First Subpop---------------------------------
	
	i = 0;	//subpop number
	popSize = p.popsize[i];
	count = 0;
	while(count < popSize)
	//for (count = 0 ; count < popSize; count++)
	{
																						#if PRINT_RECOMB
																						printf("!!!!!NEW PICK!!!!!\n");
																						#endif
		rInd1 = (int) (popSize * ran1(&gseed));
		rInd2 = (int) (popSize * ran1(&gseed));
				
																							#if PRINT_SELECTION
																							printf("\nsubpop %d -----------------------------------------------------%d\n",i,  count);
																							
																							printf("Ind1(parent : %5d) :", rInd1);
																							for (j=0;j<p.nSeq_block;j++){
																								PrintBinary(parentPop[i].indiv[rInd1].seq[j]);
																								printf(" ");}
																							printf("\nInd2(parent : %5d) :", rInd2);
																							for (j=0;j<p.nSeq_block;j++){
																								PrintBinary(parentPop[i].indiv[rInd2].seq[j]);
																								printf(" ");}
																							printf("\n");
																						//	printf("> %lf		%lf \n	", PhenoOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, selOpt[i][gen_period], 1,p), PhenoOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, selOpt[i][gen_period], 1,p));
																							
																							#endif
		#if noStorageMode
		fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
		fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
			
		#else
		fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
		fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
		#endif
		
					
					
		//	printf("++++++++++gen %d / before recom / %dth indiv+++++++++++++\n", gen, count);	
		//Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), p);
		//printf("++++++++++gen %d / after recom / %dth indiv+++++++++++++\n", gen, count);
			

		if(ran1(&gseed) < fit1 )
		//	ran1 < relative fitness
		// = ran1 * mostFit < indFit
		// = ran1 * 1 < indFit  (since mostFit = 1, see GetFitnessFromPheno)
		{
			 
			
																				#if PRINT_RECOMB
																				printf("\nBefore Recombinant (%p / %p) : ", &recombinant, recombinant.seq);
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf("\n");
																				#endif
			
			
			//Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), p);			
			Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), 0,  &recombinant, p);	
																				#if PRINT_RECOMB
																				printf("After  Recombinant (%p / %p) : ", &recombinant, recombinant.seq);
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf("\n");
																				#endif
																				
			
			
			//offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[rInd1], p);
			offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], recombinant, p);
			countHaplotype(haploTable, offspringPop[i].indiv[count].seq, p);
			
			
																					#if PRINT_SELECTION
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(parentPop[i].indiv[rInd1].seq[j]);
																					printf(" ");
																				}
																				printf(" parantIn1 \n");
																				
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf(" recombina \n---------------------------\n");
																				
																				
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(offspringPop[i].indiv[count].seq[j]);
																					printf(" ");
																				}
																				printf(" offspring \n===========================\n");
																				
																				
																				#endif 
																																			
																				
			count += 1;
		
			if (count == popSize) 
			{
																					#if PRINT_SELECTION
																					printf("\noffspringPop\n");
																					PrintPop(offspringPop, p);
																					#endif
				break;			 
			}
																					
		}
		
				
		if(ran1(&gseed) < fit2 )
		//	ran1 < relative fitness
		// = ran1 * mostFit < indFit
		// = ran1 * 1 < indFit  (since mostFit = 1, see GetFitnessFromPheno)
		{
			 
		


																				#if PRINT_RECOMB
																				printf("\nBefore Recombinant (%p / %p) : ", &recombinant, recombinant.seq);
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf("\n");

																				#endif
																				
			
			//Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), p);			
			Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), 1,  &recombinant, p);	
																				#if PRINT_RECOMB
																				printf("After  Recombinant (%p / %p) : ", &recombinant, recombinant.seq);
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf("\n");
																				#endif


			
			//offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[rInd2], p);
			offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], recombinant, p);
			countHaplotype(haploTable, offspringPop[i].indiv[count].seq, p);
			
																				#if PRINT_SELECTION
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(parentPop[i].indiv[rInd2].seq[j]);
																					printf(" ");
																				}
																				printf(" parantIn2 \n");
																				
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf(" recombina \n---------------------------\n");
																				
																				
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(offspringPop[i].indiv[count].seq[j]);
																					printf(" ");
																				}
																				printf(" offspring \n===========================\n");
																				
																				
																				#endif 
																				
																				
			count += 1;
																					
																					
			if (count == popSize) 
			{
																					#if PRINT_SELECTION
																					printf("\noffspringPop\n");
																					PrintPop(offspringPop, p);
																					#endif
				break;			 
			}
																					
		}
		
	}				
				
				
	//Second Subpop---------------------------------
	//just rearrange the order of individuals
	

	i = 1;	//subpop number
	popSize = p.popsize[i];
	
	ranArray = (int*) malloc (popSize * sizeof(int) );		//for random sorting of second subpopulation 
	for ( j = 0; j < popSize; j++)
	{
		ranArray[j] = j;
	}
	
	

	for (count = 0 ; count < popSize; count++)
	{
		//randomly arrange index number 
		//randInd = ran1(&gseed) * (popSize - count) + count ; 
		randInd = (int) (ran1(&gseed) * (popSize - count)) + count ;
			

		tmp = ranArray[count];
		ranArray[count] = ranArray[randInd];
		ranArray[randInd] = tmp;
		
		
																				#if PRINT_SELECTION
																				printf("%d -> %d | ", randInd, ranArray[count]);
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(parentPop[i].indiv[ranArray[count]].seq[j]);
																					printf(" ");
																				}
																				printf("    sub%d  indiv %d\n", i, ranArray[count]);
																				#endif
	
		
																				
																				
		
		
		
		
		
		offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[ranArray[count]], p);
		countHaplotype(haploTable, offspringPop[i].indiv[count].seq, p);
		
		//count +=1;
		
		if (count == popSize) 
		{
			break;			 
		}
			
	}
	
	

			
	free(ranArray);
	free(recombinant.seq);
	
	
	return offspringPop;
}


struct Population* SelectionAndReproducation_Guli_Refuge(struct Population* currentPop, struct Population* pop1, struct Population* pop2,unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen, list* haploTable, para_t p)
//written 171016 -> refuge
{
	/* 	Based on Wright-Fisher model
		
		In subpopulation 1, 
		  randomly select two individuals and get fitness
		  if ranNum > fitness, then the individual undergoes recombination with the other indiv - for both indiv
		  offspring will be recombinant of the individual and haploType is counted
		  repeat until count == popSize ( the population size of subPop1)
		
		In subpopulation 2 (SeedBank), 
		  randomly select two individuals 
		  individual undergoes recombination with the other indiv - for both indiv
		  offspring will be recombinant of the inidicual and haploType is counted
		  repeat until count == popSize (the population size of subPop1)
		  => everything is same with subpop1, except that this lacks checking ranNum>fit step
	*/	  
	
	
	#if PRINT_WHERE
	printf("\n\n==========SelectionAndReproducation==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========SelectionAndReproducation==========\n");
	#endif
	
	
	struct Population* parentPop;
	struct Population* offspringPop;
	int rInd1, rInd2, randInd;
	int i,j, count, k, tmp;
	int gen_period = gen % (p.period);
	double fit1, fit2;
	double oldpheno1=-100000,oldpheno2=-100000,  nowpheno;
	int phenoCount1=0, phenoCount2=0;
	int* seq;
	int* ranArray;
	int popSize;

	
	if (currentPop == pop1)								//pop1 = parentPop, pop2 = offspringPop
	{
		parentPop = pop1;
		offspringPop = pop2;
	}
	else if (currentPop == pop2)						//pop2 = parentPop, pop1 = offspringPop	
	{
		parentPop = pop2;
		offspringPop = pop1;		
	}
	
	
	struct Individual recombinant;
	recombinant.seq = (unsigned int*) malloc (p.nSeq_block * sizeof(unsigned int));
	
	
	
	//First Subpop---------------------------------
	
	i = 0;	//subpop number
	popSize = p.popsize[i];
	count = 0;
	while(count < popSize)
	//for (count = 0 ; count < popSize; count++)
	{
																						#if PRINT_RECOMB
																						printf("!!!!!NEW PICK!!!!!\n");
																						#endif
		rInd1 = (int) (popSize * ran1(&gseed));
		rInd2 = (int) (popSize * ran1(&gseed));
				
																							#if PRINT_SELECTION
																							printf("\nsubpop %d -----------------------------------------------------%d\n",i,  count);
																							
																							printf("Ind1(parent : %5d) :", rInd1);
																							for (j=0;j<p.nSeq_block;j++){
																								PrintBinary(parentPop[i].indiv[rInd1].seq[j]);
																								printf(" ");}
																							printf("\nInd2(parent : %5d) :", rInd2);
																							for (j=0;j<p.nSeq_block;j++){
																								PrintBinary(parentPop[i].indiv[rInd2].seq[j]);
																								printf(" ");}
																							printf("\n");
																						//	printf("> %lf		%lf \n	", PhenoOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, selOpt[i][gen_period], 1,p), PhenoOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, selOpt[i][gen_period], 1,p));
																							
																							#endif
		#if noStorageMode
		fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
		fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
			
		#else
		fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
		fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
		#endif
		
					
					
		//	printf("++++++++++gen %d / before recom / %dth indiv+++++++++++++\n", gen, count);	
		//Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), p);
		//printf("++++++++++gen %d / after recom / %dth indiv+++++++++++++\n", gen, count);
			

		if(ran1(&gseed) < fit1 )
		//	ran1 < relative fitness
		// = ran1 * mostFit < indFit
		// = ran1 * 1 < indFit  (since mostFit = 1, see GetFitnessFromPheno)
		{
																				#if PRINT_RECOMB
																				printf("\nBefore Recombinant (%p / %p) : ", &recombinant, recombinant.seq);for (j = 0; j < p.nSeq_block; j++){PrintBinary(recombinant.seq[j]);printf(" ");}printf("\n");
																				#endif
				
			Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), 0,  &recombinant, p);	
																				#if PRINT_RECOMB
																				printf("After  Recombinant (%p / %p) : ", &recombinant, recombinant.seq);for (j = 0; j < p.nSeq_block; j++){PrintBinary(recombinant.seq[j]);printf(" ");}printf("\n");
																				#endif
																				
			offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], recombinant, p);
			countHaplotype(haploTable, offspringPop[i].indiv[count].seq, p);
			
			
																					#if PRINT_SELECTION
																					for (j = 0; j < p.nSeq_block; j++){PrintBinary(parentPop[i].indiv[rInd1].seq[j]);printf(" ");}printf(" parantIn1 \n");
																					for (j = 0; j < p.nSeq_block; j++){PrintBinary(recombinant.seq[j]);printf(" ");}printf(" recombina \n---------------------------\n");
																					for (j = 0; j < p.nSeq_block; j++){PrintBinary(offspringPop[i].indiv[count].seq[j]);printf(" ");}printf(" offspring \n===========================\n");
																					#endif 
			count += 1;
		
			if (count == popSize) 
			{
																					#if PRINT_SELECTION
																					printf("\noffspringPop\n");PrintPop(offspringPop, p);
																					#endif
				break;			 
			}
																					
		}
		
				
		if(ran1(&gseed) < fit2 )
		//	ran1 < relative fitness
		// = ran1 * mostFit < indFit
		// = ran1 * 1 < indFit  (since mostFit = 1, see GetFitnessFromPheno)
		{
																				#if PRINT_RECOMB
																				printf("\nBefore Recombinant (%p / %p) : ", &recombinant, recombinant.seq);for (j = 0; j < p.nSeq_block; j++){PrintBinary(recombinant.seq[j]);printf(" ");}printf("\n");
																				#endif
			Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), 1,  &recombinant, p);	
																				#if PRINT_RECOMB
																				printf("After  Recombinant (%p / %p) : ", &recombinant, recombinant.seq);for (j = 0; j < p.nSeq_block; j++){PrintBinary(recombinant.seq[j]);printf(" ");}printf("\n");
																				#endif
			offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], recombinant, p);
			countHaplotype(haploTable, offspringPop[i].indiv[count].seq, p);
																					#if PRINT_SELECTION
																					for (j = 0; j < p.nSeq_block; j++){PrintBinary(parentPop[i].indiv[rInd1].seq[j]);printf(" ");}printf(" parantIn1 \n");
																					for (j = 0; j < p.nSeq_block; j++){PrintBinary(recombinant.seq[j]);printf(" ");}printf(" recombina \n---------------------------\n");
																					for (j = 0; j < p.nSeq_block; j++){PrintBinary(offspringPop[i].indiv[count].seq[j]);printf(" ");}printf(" offspring \n===========================\n");
																					#endif 
			count += 1;
																					
			if (count == popSize) 
			{
																					#if PRINT_SELECTION
																					printf("\noffspringPop\n");PrintPop(offspringPop, p);
																					#endif
				break;			 
			}
																					
		}
		
	}				
				
				
	//Second Subpop---------------------------------
	//not check fitness
	i = 1;	//subpop number
	popSize = p.popsize[i];
	count = 0;
	
	while(count < popSize)
	{
																						#if PRINT_RECOMB
																						printf("!!!!!NEW PICK!!!!!\n");
																						#endif
		rInd1 = (int) (popSize * ran1(&gseed));
		rInd2 = (int) (popSize * ran1(&gseed));
		
		
		// for rInd1----------------------	
																				#if PRINT_RECOMB
																				printf("\nBefore Recombinant (%p / %p) : ", &recombinant, recombinant.seq);for (j = 0; j < p.nSeq_block; j++){PrintBinary(recombinant.seq[j]);printf(" ");}printf("\n");
																				#endif
																				
																				
		Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), 0,  &recombinant, p);	
																				#if PRINT_RECOMB
																				printf("After  Recombinant (%p / %p) : ", &recombinant, recombinant.seq);for (j = 0; j < p.nSeq_block; j++){PrintBinary(recombinant.seq[j]);printf(" ");}printf("\n");
																				#endif

		offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], recombinant, p);
		countHaplotype(haploTable, offspringPop[i].indiv[count].seq, p);
			
																					#if PRINT_SELECTION
																					for (j = 0; j < p.nSeq_block; j++){PrintBinary(parentPop[i].indiv[rInd1].seq[j]);printf(" ");}printf(" parantIn1 \n");
																					for (j = 0; j < p.nSeq_block; j++){PrintBinary(recombinant.seq[j]);printf(" ");}printf(" recombina \n---------------------------\n");
																					for (j = 0; j < p.nSeq_block; j++){PrintBinary(offspringPop[i].indiv[count].seq[j]);printf(" ");}printf(" offspring \n===========================\n");
																					#endif 
																		
		count += 1;
		
		if (count == popSize) 
		{
																					#if PRINT_SELECTION
																					printf("\noffspringPop\n");PrintPop(offspringPop, p);
																					#endif
			break;			 
		}
		
		
		// for rInd2----------------------	
																				#if PRINT_RECOMB
																				printf("\nBefore Recombinant (%p / %p) : ", &recombinant, recombinant.seq);for (j = 0; j < p.nSeq_block; j++){PrintBinary(recombinant.seq[j]);printf(" ");}printf("\n");
																				#endif
			Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), 1,  &recombinant, p);	
																				#if PRINT_RECOMB
																				printf("After  Recombinant (%p / %p) : ", &recombinant, recombinant.seq);for (j = 0; j < p.nSeq_block; j++){PrintBinary(recombinant.seq[j]);printf(" ");}printf("\n");
																				#endif
			offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], recombinant, p);
			countHaplotype(haploTable, offspringPop[i].indiv[count].seq, p);
																					#if PRINT_SELECTION
																					for (j = 0; j < p.nSeq_block; j++){PrintBinary(parentPop[i].indiv[rInd1].seq[j]);printf(" ");}printf(" parantIn1 \n");
																					for (j = 0; j < p.nSeq_block; j++){PrintBinary(recombinant.seq[j]);printf(" ");}printf(" recombina \n---------------------------\n");
																					for (j = 0; j < p.nSeq_block; j++){PrintBinary(offspringPop[i].indiv[count].seq[j]);printf(" ");}printf(" offspring \n===========================\n");
																					#endif 
			count += 1;
																					
			if (count == popSize) 
			{
																					#if PRINT_SELECTION
																					printf("\noffspringPop\n");PrintPop(offspringPop, p);
																					#endif
				break;			 
			}
	}			
	
	

			

	free(recombinant.seq);
	
	
	return offspringPop;
}

struct Population* SelectionAndReproducation_Guli_bothSelect(struct Population* currentPop, struct Population* pop1, struct Population* pop2,unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen, list* haploTable, para_t p)
{
	/* 	Based on Wright-Fisher model
		
		In subpopulation 1, 
		  randomly select two individuals and get fitness
		  if ranNum > fitness, then the individual undergoes recombination with the other indiv - for both indiv
		  offspring will be recombinant of the individual and haploType is counted
		  repeat until count == popSize ( the population size of subPop1)
		
		In subpopulation 2 (SeedBank), 
		  the order of individual is randomly reassigned in offspring population using for loop
		  while copying indivs to offpop, haploType is counted
	*/	  
	
	
	#if PRINT_WHERE
	printf("\n\n==========SelectionAndReproducation==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========SelectionAndReproducation==========\n");
	#endif
	
	
	struct Population* parentPop;
	struct Population* offspringPop;
	int rInd1, rInd2, randInd;
	int i,j, count, k, tmp;
	int gen_period = gen % (p.period);
	double fit1, fit2;
	double oldpheno1=-100000,oldpheno2=-100000,  nowpheno;
	int phenoCount1=0, phenoCount2=0;
	int* seq;
	int* ranArray;
	int popSize;
	
	
	if (currentPop == pop1)								//pop1 = parentPop, pop2 = offspringPop
	{
		parentPop = pop1;
		offspringPop = pop2;
	}
	else if (currentPop == pop2)						//pop2 = parentPop, pop1 = offspringPop	
	{
		parentPop = pop2;
		offspringPop = pop1;		
	}

	struct Individual recombinant;
	recombinant.seq = (unsigned int*) malloc (p.nSeq_block * sizeof(unsigned int));
	
	
	
	//First Subpop---------------------------------
	
	i = 0;	//subpop number
	popSize = p.popsize[i];
	count = 0;
	while(count < popSize)
	//for (count = 0 ; count < popSize; count++)
	{
																						#if PRINT_RECOMB
																						printf("!!!!!NEW PICK!!!!!\n");
																						#endif
		rInd1 = (int) (popSize * ran1(&gseed));
		rInd2 = (int) (popSize * ran1(&gseed));
				
																							#if PRINT_SELECTION
																							printf("\nsubpop %d -----------------------------------------------------%d\n",i,  count);
																							
																							printf("Ind1(parent : %5d) :", rInd1);
																							for (j=0;j<p.nSeq_block;j++){
																								PrintBinary(parentPop[i].indiv[rInd1].seq[j]);
																								printf(" ");}
																							printf("\nInd2(parent : %5d) :", rInd2);
																							for (j=0;j<p.nSeq_block;j++){
																								PrintBinary(parentPop[i].indiv[rInd2].seq[j]);
																								printf(" ");}
																							printf("\n");
																						//	printf("> %lf		%lf \n	", PhenoOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, selOpt[i][gen_period], 1,p), PhenoOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, selOpt[i][gen_period], 1,p));
																							
																							#endif
		#if noStorageMode
		fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
		fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
			
		#else
		fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
		fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
		#endif
		
					
					
		//	printf("++++++++++gen %d / before recom / %dth indiv+++++++++++++\n", gen, count);	
		//Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), p);
		//printf("++++++++++gen %d / after recom / %dth indiv+++++++++++++\n", gen, count);
			

		if(ran1(&gseed) < fit1 )
		//	ran1 < relative fitness
		// = ran1 * mostFit < indFit
		// = ran1 * 1 < indFit  (since mostFit = 1, see GetFitnessFromPheno)
		{
			 
			
																				#if PRINT_RECOMB
																				printf("\nBefore Recombinant (%p / %p) : ", &recombinant, recombinant.seq);
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf("\n");
																				#endif
			
			
			//Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), p);			
			Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), 0,  &recombinant, p);	
																				#if PRINT_RECOMB
																				printf("After  Recombinant (%p / %p) : ", &recombinant, recombinant.seq);
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf("\n");
																				#endif
																				
			
			
			//offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[rInd1], p);
			offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], recombinant, p);
			countHaplotype(haploTable, offspringPop[i].indiv[count].seq, p);
			
			
																					#if PRINT_SELECTION
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(parentPop[i].indiv[rInd1].seq[j]);
																					printf(" ");
																				}
																				printf(" parantIn1 \n");
																				
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf(" recombina \n---------------------------\n");
																				
																				
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(offspringPop[i].indiv[count].seq[j]);
																					printf(" ");
																				}
																				printf(" offspring \n===========================\n");
																				
																				
																				#endif 
																																			
																				
			count += 1;
		
			if (count == popSize) 
			{
																					#if PRINT_SELECTION
																					printf("\noffspringPop\n");
																					PrintPop(offspringPop, p);
																					#endif
				break;			 
			}
																					
		}
		
				
		if(ran1(&gseed) < fit2 )
		//	ran1 < relative fitness
		// = ran1 * mostFit < indFit
		// = ran1 * 1 < indFit  (since mostFit = 1, see GetFitnessFromPheno)
		{
			 
		


																				#if PRINT_RECOMB
																				printf("\nBefore Recombinant (%p / %p) : ", &recombinant, recombinant.seq);
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf("\n");

																				#endif
																				
			
			//Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), p);			
			Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), 1,  &recombinant, p);	
																				#if PRINT_RECOMB
																				printf("After  Recombinant (%p / %p) : ", &recombinant, recombinant.seq);
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf("\n");
																				#endif


			
			//offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[rInd2], p);
			offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], recombinant, p);
			countHaplotype(haploTable, offspringPop[i].indiv[count].seq, p);
			
																				#if PRINT_SELECTION
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(parentPop[i].indiv[rInd2].seq[j]);
																					printf(" ");
																				}
																				printf(" parantIn2 \n");
																				
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf(" recombina \n---------------------------\n");
																				
																				
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(offspringPop[i].indiv[count].seq[j]);
																					printf(" ");
																				}
																				printf(" offspring \n===========================\n");
																				
																				
																				#endif 
																				
																				
			count += 1;
																					
																					
			if (count == popSize) 
			{
																					#if PRINT_SELECTION
																					printf("\noffspringPop\n");
																					PrintPop(offspringPop, p);
																					#endif
				break;			 
			}
																					
		}
		
	}				
				
				
	//Second Subpop---------------------------------
	
	i = 1;	//subpop number
	popSize = p.popsize[i];
	count = 0;
	while(count < popSize)
	//for (count = 0 ; count < popSize; count++)
	{
																						#if PRINT_RECOMB
																						printf("!!!!!NEW PICK!!!!!\n");
																						#endif
		rInd1 = (int) (popSize * ran1(&gseed));
		rInd2 = (int) (popSize * ran1(&gseed));
				
																							#if PRINT_SELECTION
																							printf("\nsubpop %d -----------------------------------------------------%d\n",i,  count);
																							
																							printf("Ind1(parent : %5d) :", rInd1);
																							for (j=0;j<p.nSeq_block;j++){
																								PrintBinary(parentPop[i].indiv[rInd1].seq[j]);
																								printf(" ");}
																							printf("\nInd2(parent : %5d) :", rInd2);
																							for (j=0;j<p.nSeq_block;j++){
																								PrintBinary(parentPop[i].indiv[rInd2].seq[j]);
																								printf(" ");}
																							printf("\n");
																						//	printf("> %lf		%lf \n	", PhenoOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, selOpt[i][gen_period], 1,p), PhenoOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, selOpt[i][gen_period], 1,p));
																							
																							#endif
		#if noStorageMode
		fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
		fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
			
		#else
		fit1 = FitOfIndiv(&(parentPop[i].indiv[rInd1]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
		fit2 = FitOfIndiv(&(parentPop[i].indiv[rInd2]), nonsyn, delta, pheno_wild, selOpt[i][gen_period],p);
		#endif
		
					
					
		//	printf("++++++++++gen %d / before recom / %dth indiv+++++++++++++\n", gen, count);	
		//Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), p);
		//printf("++++++++++gen %d / after recom / %dth indiv+++++++++++++\n", gen, count);
			

		if(ran1(&gseed) < fit1 )
		//	ran1 < relative fitness
		// = ran1 * mostFit < indFit
		// = ran1 * 1 < indFit  (since mostFit = 1, see GetFitnessFromPheno)
		{
			 
			
																				#if PRINT_RECOMB
																				printf("\nBefore Recombinant (%p / %p) : ", &recombinant, recombinant.seq);
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf("\n");
																				#endif
			
			
			//Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), p);			
			Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), 0,  &recombinant, p);	
																				#if PRINT_RECOMB
																				printf("After  Recombinant (%p / %p) : ", &recombinant, recombinant.seq);
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf("\n");
																				#endif
																				
			
			
			//offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[rInd1], p);
			offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], recombinant, p);
			countHaplotype(haploTable, offspringPop[i].indiv[count].seq, p);
			
			
																					#if PRINT_SELECTION
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(parentPop[i].indiv[rInd1].seq[j]);
																					printf(" ");
																				}
																				printf(" parantIn1 \n");
																				
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf(" recombina \n---------------------------\n");
																				
																				
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(offspringPop[i].indiv[count].seq[j]);
																					printf(" ");
																				}
																				printf(" offspring \n===========================\n");
																				
																				
																				#endif 
																																			
																				
			count += 1;
		
			if (count == popSize) 
			{
																					#if PRINT_SELECTION
																					printf("\noffspringPop\n");
																					PrintPop(offspringPop, p);
																					#endif
				break;			 
			}
																					
		}
		
				
		if(ran1(&gseed) < fit2 )
		//	ran1 < relative fitness
		// = ran1 * mostFit < indFit
		// = ran1 * 1 < indFit  (since mostFit = 1, see GetFitnessFromPheno)
		{
			 
		


																				#if PRINT_RECOMB
																				printf("\nBefore Recombinant (%p / %p) : ", &recombinant, recombinant.seq);
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf("\n");

																				#endif
																				
			
			//Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), p);			
			Recombination(&(parentPop[i].indiv[rInd1]), &(parentPop[i].indiv[rInd2]), 1,  &recombinant, p);	
																				#if PRINT_RECOMB
																				printf("After  Recombinant (%p / %p) : ", &recombinant, recombinant.seq);
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf("\n");
																				#endif


			
			//offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], parentPop[i].indiv[rInd2], p);
			offspringPop[i].indiv[count] = CopyIndividual(offspringPop[i].indiv[count], recombinant, p);
			countHaplotype(haploTable, offspringPop[i].indiv[count].seq, p);
			
																				#if PRINT_SELECTION
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(parentPop[i].indiv[rInd2].seq[j]);
																					printf(" ");
																				}
																				printf(" parantIn2 \n");
																				
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(recombinant.seq[j]);
																					printf(" ");
																				}
																				printf(" recombina \n---------------------------\n");
																				
																				
																				for (j = 0; j < p.nSeq_block; j++)
																				{
																					PrintBinary(offspringPop[i].indiv[count].seq[j]);
																					printf(" ");
																				}
																				printf(" offspring \n===========================\n");
																				
																				
																				#endif 
																				
																				
			count += 1;
																					
																					
			if (count == popSize) 
			{
																					#if PRINT_SELECTION
																					printf("\noffspringPop\n");
																					PrintPop(offspringPop, p);
																					#endif
				break;			 
			}
																					
		}
		
	}				
		


	free(recombinant.seq);
	
	
	return offspringPop;
}

void Mutation(struct Population* pop, para_t p)
//1114 version
{
	//mutation occurs regardless of the syn / nonsyn
	#if PRINT_WHERE
	printf("\n\n==========Mutation==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========Mutation==========\n");
	#endif
	
	int nMut;
	int rNum, rSite;
	unsigned int seg, binary=0;
	int s, m;
	struct Individual* ranIdv;
	
	

	nMut = poidev(p.Nu*_block_*p.nSeq_block, &gseed);
	
	
	for ( m = 0; m < nMut ; m++)
	{
		#ifdef _MUT12_ 
		//mutation occurs in both subpopulation
		
		rNum = (int) (ran1(&gseed) * p.nIndT);
		if (rNum < p.popsize[0])
			ranIdv = &(pop[0].indiv[rNum]);
		else
			ranIdv = &(pop[1].indiv[rNum-p.popsize[0]]);
	
		
		
		#elif defined(_MUT1_)
		//mutation occurs in only sub0
		rNum = (int) (ran1(&gseed) * p.popsize[0]);
		ranIdv = &(pop[0].indiv[rNum]);
		
		
		#else
		printf("++++++++++++++++++++++++++++++\n+[Warning] Please compile with option+\n++++++++++++++++++++++++++++++\n");
		#endif
		
		
		rSite = (int) (ran1(&gseed)*p.nSeq);	
		
											#if PRINT_SEQ
											if (rNum < p.popsize[0])
												printf(" In sub0, indiv %d" ,rNum);
											else	
												printf(" In sub1, indiv %d" ,rNum-p.popsize[0]);
											printf(" 	block %d site %d\n", rSite/_block_, rSite%_block_);
											#endif
		
		
		seg = ranIdv->seq[rSite/_block_];
											#if PRINT_SEQ
											PrintBinary(seg);
											printf("\n");
											#endif
		
		seg = seg ^ (1 << ((rSite%_block_)));
											#if PRINT_SEQ
											PrintBinary_asterik(1 << ((rSite%_block_)));
											printf("\n");
											
											PrintBinary(seg);
											printf("\n\n");
											#endif
		
		ranIdv->seq[rSite/_block_] = seg;
											
											#if PRINT_SEQ
											PrintPop(pop, p);
											#endif
	}
}




/* mutation last version
void Mutation(struct Population* pop, para_t p)
{
	//mutation occurs regardless of the syn / nonsyn
	#if PRINT_WHERE
	printf("\n\n==========Mutation==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========Mutation==========\n");
	#endif
	
	int nMut;
	int rNum, rSite;
	unsigned int seg, binary=0;
	int s, m;
	struct Individual* ranIdv;
	
	

	nMut = poidev(p.Nu*_block_*p.nSeq_block*2, &gseed);
	for ( m = 0; m < nMut ; m++)
	{
		rNum = (int) (ran1(&gseed) * p.nIndT);
		if (rNum < p.popsize[0])
			ranIdv = &(pop[0].indiv[rNum]);
		else
			ranIdv = &(pop[1].indiv[rNum-p.popsize[0]]);
			
		rSite = (int) (ran1(&gseed)*p.nSeq);	
		
											#if PRINT_SEQ
											if (rNum < p.popsize[0])
												printf(" In sub0, indiv %d" ,rNum);
											else	
												printf(" In sub1, indiv %d" ,rNum-p.popsize[0]);
											printf(" 	block %d site %d\n", rSite/_block_, rSite%_block_);
											#endif
		
		
		seg = ranIdv->seq[rSite/_block_];
											#if PRINT_SEQ
											PrintBinary(seg);
											printf("\n");
											#endif
		
		seg = seg ^ (1 << ((rSite%_block_)));
											#if PRINT_SEQ
											PrintBinary_asterik(1 << ((rSite%_block_)));
											printf("\n");
											
											PrintBinary(seg);
											printf("\n\n");
											#endif
		
		ranIdv->seq[rSite/_block_] = seg;
											
											#if PRINT_SEQ
											PrintPop(pop, p);
											#endif
	}
}

void Mutation_SeedBank(struct Population* pop, para_t p)
{
	//mutation occurs regardless of the syn / nonsyn
	#if PRINT_WHERE
	printf("\n\n==========Mutation==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========Mutation==========\n");
	#endif
	
	int nMut;
	int rNum, rSite;
	unsigned int seg, binary=0;
	int s, m;
	struct Individual* ranIdv;
											
										
	

	nMut = poidev(p.Nu*_block_*p.nSeq_block, &gseed);
											
											#if PRINT_SEQ
											printf(" < nMut = %d >\n", nMut);
											#endif
	for ( m = 0; m < nMut ; m++)
	{
		rNum = (int) (ran1(&gseed) * p.popsize[0]);
		ranIdv = &(pop[0].indiv[rNum]);
		
		
		rSite = (int) (ran1(&gseed)*p.nSeq);	
		
											#if PRINT_SEQ
											printf(" In sub0, indiv %d" ,rNum);
											printf(" 	block %d site %d\n", rSite/_block_, rSite%_block_);
											#endif
		
		
		seg = ranIdv->seq[rSite/_block_];
											#if PRINT_SEQ
											PrintBinary(seg);
											printf("\n");
											#endif
		
		seg = seg ^ (1 << ((rSite%_block_)));
											#if PRINT_SEQ
											PrintBinary_asterik(1 << ((rSite%_block_)));
											printf("\n");
											
											PrintBinary(seg);
											printf("\n\n");
											#endif
		
		ranIdv->seq[rSite/_block_] = seg;
											
											#if PRINT_SEQ
											PrintPop(pop, p);
											#endif
	}
}
*/





void PrintFreq(struct Population* pop, int gen, para_t p)
{
	//was GetHet originally
	
	#if PRINT_WHERE
	printf("\n\n==========PrintFreq==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========PrintFreq==========\n");
	#endif
	
	
	int i, j, k, check = 1;
	char bit;
	unsigned int input;
	int countAllele[_block_];
	fprintf(g_freqfile,"%d ", gen);
	
	
	for (j = 0; j < p.nSeq_block; j++)
	{
		memset(countAllele, 0, sizeof(int)*_block_);
		for (k = 0; k < p.popsize[0]; k++)				//subpop0
		{
			input = pop[0].indiv[k].seq[j];
			
			for (i = _block_-1; i>=0 ; i--)
			{
				bit = ( input & ( 1 << i ))?1:0;
				if (bit ==1)
					countAllele[i]++;
			}
		}
		for (k = 0; k < p.popsize[1]; k++)				//subpop1
		{
			input = pop[1].indiv[k].seq[j];
			
			for (i = _block_-1; i>=0 ; i--)
			{
				bit = ( input & ( 1 << i ))?1:0;
				if (bit ==1)
					countAllele[i]++;
			}
		}
	
		
	
		for (i = _block_-1; i>=0 ; i--)
		{
			fprintf(g_freqfile, "%d ", countAllele[i]);
		}
		fprintf(g_freqfile," /");
		
	}
	fprintf(g_freqfile,"\n");
	
	
}



int PrintFreqAndFixLossCheck(struct Population* pop, int gen, para_t p)
{
	//was GetHet originally
	
	#if PRINT_WHERE
	printf("\n\n==========PrintFreq==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========PrintFreq==========\n");
	#endif
	
	
	
	
	int i, j, k, check = 1;
	char bit;
	unsigned int input;
	int countAllele[_block_];
	fprintf(g_freqfile,"%d ", gen);
	
	int selecSeg, selecSite;
	selecSeg = g_rSite/_block_;
	selecSite = g_rSite%_block_;
	int lossOrfix;
	
	for (j = 0; j < p.nSeq_block; j++)
	{
		memset(countAllele, 0, sizeof(int)*_block_);
		for (k = 0; k < p.popsize[0]; k++)				//subpop0
		{
			input = pop[0].indiv[k].seq[j];
			
			for (i = _block_-1; i>=0 ; i--)
			{
				bit = ( input & ( 1 << i ))?1:0;
				if (bit ==1)
					countAllele[i]++;
			}
		}
		for (k = 0; k < p.popsize[1]; k++)				//subpop1
		{
			input = pop[1].indiv[k].seq[j];
			
			for (i = _block_-1; i>=0 ; i--)
			{
				bit = ( input & ( 1 << i ))?1:0;
				if (bit ==1)
					countAllele[i]++;
			}
		}
	
		
	
		for (i = _block_-1; i>=0 ; i--)
		{
			fprintf(g_freqfile, "%d ", countAllele[i]);
		}
		fprintf(g_freqfile," /");
		
		
		//check lossOrfix
		if (j == selecSeg)
		{
			if (countAllele[selecSite] < 1)
				lossOrfix = 1;
			else if ( countAllele[selecSite] > p.nIndT - 1)
				lossOrfix = 1;
			
			else
				lossOrfix = 0;
		}
			
		//END check lossOrfix
			
		
	}
	fprintf(g_freqfile,"\n");
	
	
	return lossOrfix;
	
	
}










void main(int argc, char *argv[])
{
	const char* inpname;
	const char* outname;
	const char* freqname;
	const char* phenoname;
	const char* haploname;
	const char* hfpname;
	char filename[100];

	para_t parameter;
		
	struct Population* pop1;
	struct Population* pop2;
	struct Population* pop;
	
	int try, gen, nGen, dm,nIndTotal;
	int periodCount;
	int i, j;
	int  nNonCountTry=0, lossOrfix, maxGen=-1, nTry_BS=0, nTry_noBS=0;
	double pheno_wild;
	
	unsigned int* template_nonsys;
	double* delta;
	double** selOpt;
	double** selOpt_origianl;
	double** freqMut;
	double* recomRateAtSite;	
	
	
	/*time.h*/
	clock_t before;
	double clockResult;
	before = clock();

	//get name of input file and fileTitle(hfpname) for .out, .freq, .pheno file
	if(argc==3){
		//thisSeed = atoi(argv[3]);
		//thisSeed = 0;
        inpname = argv[1];
		hfpname = argv[2];
		
   
    }
	else{
		printf("./PROGRAM_NAME inpname hfpname");
	}	
	sprintf(filename, "%s.out",hfpname);
	outname = filename;
	


	//Init_Simulation
	parameter = ReadParameter(inpname, outname, parameter);
	gseed = (-1) * parameter.seed;
	nGen = parameter.nGen;
	
	
	
	//GetArray
	template_nonsys = (unsigned int*) malloc ( parameter.nSeq_block * sizeof(unsigned int));
	delta = (double*) malloc (sizeof(double) * parameter.nNon);
	
	selOpt = (double**) malloc (sizeof(double*) * parameter.nDeme);
	for (i =0;i<parameter.nDeme;i++)
	{
		selOpt[i] = (double*) malloc (sizeof(double) * parameter.period) ;				
	}
	
	selOpt_origianl = (double**) malloc (sizeof(double*) * parameter.nDeme);
	for (i =0;i<parameter.nDeme;i++)
	{
		selOpt_origianl[i] = (double*) malloc (sizeof(double) * parameter.period) ;				
	}

	//GetArray End

	
	
	
	//getArray for struct Population	
	pop1 = allocPop(parameter);	
	pop2 = allocPop(parameter);
	

	//create haploTable
	list* haploTable = (list*)malloc(sizeof(list));
	haploTable->count = 0;
	haploTable->head = NULL;
	
	
	////Try
	try = 0;
	while ( (try < parameter.nTry)  && (nNonCountTry<10000) )
	{
		printf("#Try %d (%d)-----------------------------------------------------\n", try, nNonCountTry);
 
	
	//	Init_Trial
		AssignDeltaEachNon( delta, parameter);
		AssignSelOptima(selOpt_origianl, 1, parameter);
		AssignNonsynSite(template_nonsys , parameter);
		//clean Population
		for (dm=0; dm<parameter.nDeme; dm++)		
		{
			for(i = 0; i<parameter.popsize[dm]; i++)
			{
				memset(pop1[dm].indiv[i].seq, 0, sizeof(unsigned int)*parameter.nSeq_block);

			}	
		}
		
		pheno_wild = ran1(&gseed) * (parameter.sMax * 2) - parameter.sMax;
	//End of Init_Trial
		#if PRINT_TML 
		printf(" pheno_wild : %lf\n", pheno_wild);
		
		#endif


		pop = pop1;
		
		//only for balancing selection
		//pop[0].indiv[0].seq[g_rSite/_block_] = pop[0].indiv[0].seq[g_rSite/_block_] |  template_nonsys[g_rSite/_block_];
		pop[0].indiv[0].seq[g_rSite/_block_] = pop[0].indiv[0].seq[g_rSite/_block_] | (1 << (g_rSite%_block_));
		
		
				 
				 
				 
	
	
											#if PRINT_RESULT
											printf("\n\n==============Ready for new generation==============\n");
											printf("--------pop --------------\n");
											PrintPop(pop, parameter);
											printf("--------pop1--------------\n");
											PrintPop(pop1, parameter);		
											printf("--------pop2--------------\n");
											PrintPop(pop2, parameter);		
											#endif		



		//start of generation
		gen = 0;
		lossOrfix = 0;
		periodCount = 0;

		
		
		//open new freq and pheno file for new trial
		sprintf(filename, "%s_try%d.freq",hfpname, try);
		freqname = filename;
		g_freqfile = fopen(freqname, "w");
		
		sprintf(filename, "%s_try%d.pheno",hfpname, try);
		phenoname = filename;
		//g_phenofile = fopen(phenoname, "w");
		
		sprintf(filename, "%s_try%d.haplo",hfpname, try);
		haploname = filename;
		g_haplofile = fopen(haploname, "w");

		
		

		//print Trial infos to .out, .pheno, .freq
		//fprintf(g_phenofile, "Try%d\t", try);
		fprintf(g_freqfile, "Try%d | ", try);
		fprintf(g_haplofile, "Try%d | ", try);
	
		 
		
		
		
		//print to g_freqfile that  where is nonsynSite
		for (i = 0 ; i < parameter.nSeq_block; i++)
		{
			PrintBinary_out(template_nonsys[i], g_freqfile);
			fprintf(g_freqfile," ");
			
			PrintBinary_out(template_nonsys[i], g_haplofile);
			fprintf(g_haplofile," ");
			
			//fprintf(g_phenofile, "wildPheno : %lf\n", pheno_wild);
		}
		fprintf(g_haplofile,"| wildPheno %lf \n", pheno_wild);
		fprintf(g_freqfile,"\n");
		
		
	
		
		
		
		while( (gen<nGen) && (lossOrfix != 1) )
		{
			//haploTable initiation
			init(haploTable);
			
			
			//disturb SelOpt fluctuation
			if ( periodCount % parameter.period == 0)
			{
				DisturbSelOpt(selOpt, selOpt_origianl, parameter.epsilon, parameter.period, parameter.nDeme);
				periodCount = 0;
			}
			
			//fprintf(g_phenofile, "Gen%d\t", gen);
											#if PRINT_WHERE
											printf(">Generation %30d\n", gen);
											
											#else
											if (gen%(5*_PRINTOUT_EVERY_) == 0)
											{
												clockResult = (double)(clock()-before)/CLOCKS_PER_SEC;
												printf(">Generation %30d / time : %lf sec\n", gen, clockResult);
											}
											#endif


			
			//mutation introduction
			Mutation(pop, parameter);			//this mutation is 1114 version
			//Mutation_SeedBank(pop, parameter);
											#if PRINT_RESULT
											printf("--------pop--------------\n");
											PrintPop(pop, parameter);
											#endif
			
		
			//migration
			Migration(pop, parameter);
											#if PRINT_RESULT
											printf("--------pop--------------\n");
											PrintPop(pop, parameter);
											#endif
										
										
			//SelectionAndReproducation	
			#if defined(_SELC_)
			//printf("SelectionAndReproducation_Guli_bothSelect\n");
			pop = SelectionAndReproducation_Guli_bothSelect(pop, pop1, pop2,template_nonsys, delta, pheno_wild, selOpt, gen, haploTable, parameter);
			#elif defined(_SELSB_)
			//printf("SelectionAndReproducation_Guli_SeedBank\n");
			pop = SelectionAndReproducation_Guli_SeedBank(pop, pop1, pop2,template_nonsys, delta, pheno_wild, selOpt, gen, haploTable, parameter);
			#elif defined(_SELRS_)
			//printf("SelectionAndReproducation_Guli_Refuge\n");
			pop = SelectionAndReproducation_Guli_Refuge(pop, pop1, pop2,template_nonsys, delta, pheno_wild, selOpt, gen, haploTable, parameter);
			#else
			printf("++++++++++++++++++++++++++++++\n+[Warning] Please compile with option+\n++++++++++++++++++++++++++++++\n");
			#endif
			
											#if PRINT_RESULT
											printf("--------pop--------------\n");
											PrintPop(pop, parameter);
											#endif								
											
			//fprintf(g_phenofile, "\n");
			
			#if multilocus
			PrintFreq(pop, gen, parameter);
			#else
			lossOrfix = PrintFreqAndFixLossCheck(pop, gen, parameter);
			#endif
			
		
			//every generation, print haploTable to outfile 
			fprintf(g_haplofile, "gen%d : %lf\n", gen, selOpt[0][gen % parameter.period]);
			print_Haplo_outfile(haploTable, template_nonsys, delta, pheno_wild, parameter);
			//END haplotype Table
			
			
			
			gen++;	
			periodCount++;
		
		}//End of one Generation		

		
		
		//get maximum generation
		if (gen > maxGen)
			maxGen = gen;
		
		//check whether BS occurred or not
		if (gen == nGen)
		{
			
			//print deltas to outfile
			fprintf(g_outfile, "Try%d\t", try);
			for (i = 0; i<parameter.nNon; i++)
			{
				fprintf(g_outfile, "%lf ", delta[i]);
			}
			fprintf(g_outfile, "\n");
			fflush(g_outfile);
			
			
			
			
			try++;
			nNonCountTry = 0;
			nTry_BS++;
			
				
			
	
			
			
			
		}
		else 
		{
			nNonCountTry++;
			nTry_noBS++;	
			
			
			fprintf(g_freqfile, "-----no complete balancing selection------------");
			//fprintf(g_phenofile, "-----no complete balancing selection------------");
			fprintf(g_haplofile, "-----no complete balancing selection------------");
			
		}
		//END check whether BS occurred or not
		
		
		fclose(g_freqfile);
		//fclose(g_phenofile);
		fclose(g_haplofile);
		//
		//try++;
		
		
		
		
	}//End of one Trial
	
	#if multilocus
	fprintf(g_outfile,"\n------multilocus simulation : cannot check balancing selection\n");
	#else
	if (nNonCountTry == 10000)
	{
		fprintf(g_outfile, "No balancing selection with this parameter -> longest generation : %d\n",maxGen );
	}
	fprintf(g_outfile,"\n trials with / without Balancing Selection = %d / %d \n", nTry_BS, nTry_noBS);
	#endif
	
	/*time.h*/
	clockResult = (double)(clock()-before)/CLOCKS_PER_SEC;
	fprintf(g_outfile,"\n totalTime = %lf sec\n", clockResult);
	
	//fclose(g_freqfile);
	fclose(g_outfile);
	//fclose(g_phenofile);
	
	free(parameter.popsize);
	
}
	



