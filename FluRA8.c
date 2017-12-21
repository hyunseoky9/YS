/*************************************************************************
/* This code is edited from FluRA5.c on November 9th, 2016
/* This code does not allow back mutation at epitope sites
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define nNonsy3 0
#define nSynon3 0
#define nNonsy4 0
#define nSynon4 0

typedef struct{
	char *AA;
	char *Sy;
	char *AA2;
	char *Sy2;
	char *AA3;
	char *Sy3;
	char *AA4;
	char *Sy4;
	double *index;
	float *Rfit;
	int numDel1;
	int numDel2;
}
vtype;

vtype **G1, **G2;

//parameters
int Ndeme, Tmax, Tburn, Tsample, DinY, Ntry, SamInt, Nsamp, Npop, dispv, dWin, numEpitope, reduced_u;
float Kmax, r_mut, sel_co, **r_mig, mostfit;
long seed;

int nNonsy, nSynon, *Nvrs;
char *beneAllele;

int bottleneckSize, doRecomb;
float sel_d;

double ***sampAnces;
char ***sampEpi;

FILE *outfile, *infile;

void Parameters(const char* inpname, const char* outname, long thisSeed)
{
	int dm, dm2;

	outfile = fopen(outname, "w" );

	if ( (infile = fopen( inpname, "r")) == NULL )
	{
		printf( "cannot open input file\n");
		exit(-1);
	}

	fscanf(infile, "%*s%d", &nNonsy);
	fscanf(infile, "%*s%d", &nSynon);
	fscanf(infile, "%*s%d", &Ndeme);
	
	r_mig = (float **)malloc(Ndeme*sizeof(float *));
	for (dm=0; dm<Ndeme; dm++)
		r_mig[dm] = (float *)malloc(Ndeme*sizeof(float));
	
	fscanf(infile, "%*s%f", &Kmax);
	fscanf(infile, "%*s%d", &DinY);
	fscanf(infile, "%*s%d", &Tburn);
	fscanf(infile, "%*s%d", &Tsample);
	fscanf(infile, "%*s%f", &r_mut);

	fscanf(infile, "%*s");
	for (dm=0; dm<Ndeme; dm++)
	{
		for (dm2=0; dm2<Ndeme; dm2++)
		{
			fscanf(infile, "%f", r_mig[dm]+dm2);
			//printf("%f\n", r_mig[dm][dm2]);
		}
	}
	fscanf(infile, "%*s%d", &SamInt);
	fscanf(infile, "%*s%d", &Nsamp);
	fscanf(infile, "%*s%d", &dWin);
	fscanf(infile, "%*s%d", &dispv);
	fscanf(infile, "%*s%d", &Ntry);
	fscanf(infile, "%*s%ld", &seed);
	if(thisSeed!=0) {
		seed = thisSeed;
	}
	fscanf(infile, "%*s%d", &bottleneckSize);
	if(bottleneckSize==0) bottleneckSize = Kmax;
	fscanf(infile, "%*s%d", &doRecomb);

	fscanf(infile, "%*s%d", &numEpitope);
	fscanf(infile, "%*s%f", &sel_co);
	fscanf(infile, "%*s%d", &reduced_u);
	fscanf(infile, "%*s%f", &sel_d);

	fprintf(outfile, "nNonsy = %d\n", nNonsy);
	fprintf(outfile, "nSynon = %d\n", nSynon);
	fprintf(outfile, "Ndeme = %d\n", Ndeme);
	fprintf(outfile, "Kmax = %f\n", Kmax);
	fprintf(outfile, "DinY = %d\n", DinY);
	fprintf(outfile, "Tburn = %d\n", Tburn);
	fprintf(outfile, "Tsample = %d\n", Tsample);
	fprintf(outfile, "r_mut = %f\n", r_mut);
	
	fprintf(outfile, "r_mig = \n");
	for (dm=0; dm<Ndeme; dm++)
	{
		for (dm2=0; dm2<Ndeme; dm2++)
			fprintf(outfile, " %f", r_mig[dm][dm2]);
		fprintf(outfile, "\n");
	}
	fprintf(outfile, "SamInt = %d\n", SamInt);
	fprintf(outfile, "Nsamp = %d\n", Nsamp);
	fprintf(outfile, "DayWin = %d\n", dWin);
	fprintf(outfile, "DispV = %d\n", dispv);
	fprintf(outfile, "Ntry = %d\n", Ntry);
	fprintf(outfile, "numEpitope = %d\n", numEpitope);
	fprintf(outfile, "sel_co = %f\n", sel_co);
	fprintf(outfile, "seed = %ld\n", seed);
	fprintf(outfile, "reduced_u = %d\n", reduced_u);

	fprintf(outfile, "bottleneckSize = %d\n", bottleneckSize);
	fprintf(outfile, "doRecomb? %d\n", doRecomb);

	fprintf(outfile, "sel_d = %f\n", sel_d);
	seed *= -1;
}

void GetArrays()
{
	int i, dm, k, j, sTep;

	Npop = (int) (2*Kmax);
	Tmax = Tburn+Tsample;

	beneAllele = (char*) malloc(numEpitope*sizeof(char));

	G1 = (vtype **)malloc(Ndeme*sizeof(vtype *));
	G2 = (vtype **)malloc(Ndeme*sizeof(vtype *));
	for (dm=0; dm<Ndeme; dm++)
	{
		G1[dm] = (vtype *)malloc(Npop*sizeof(vtype));
		G2[dm] = (vtype *)malloc(Npop*sizeof(vtype));
		for (k=0; k<Npop; k++)
		{
			G1[dm][k].AA = (char *)malloc(nNonsy*sizeof(char));
			G1[dm][k].Sy = (char *)malloc(nSynon*sizeof(char));
			G1[dm][k].Rfit = (float *)malloc(Ndeme*sizeof(float));
			G1[dm][k].index = (double *)malloc(numEpitope*sizeof(double));
			G1[dm][k].AA2 = (char *)malloc(nNonsy*sizeof(char));
			G1[dm][k].Sy2 = (char *)malloc(nSynon*sizeof(char));
			G1[dm][k].AA3 = (char *)malloc(nNonsy3*sizeof(char));
			G1[dm][k].Sy3 = (char *)malloc(nSynon3*sizeof(char));
			G1[dm][k].AA4 = (char *)malloc(nNonsy4*sizeof(char));
			G1[dm][k].Sy4 = (char *)malloc(nSynon4*sizeof(char));


			G2[dm][k].AA =(char *)malloc(nNonsy*sizeof(char));
			G2[dm][k].Sy = (char *)malloc(nSynon*sizeof(char));
			G2[dm][k].Rfit = (float *)malloc(Ndeme*sizeof(float));
			G2[dm][k].index = (double *)malloc(numEpitope*sizeof(double));
			G2[dm][k].AA2 = (char *)malloc(nNonsy*sizeof(char));
			G2[dm][k].Sy2 = (char *)malloc(nSynon*sizeof(char));
			G2[dm][k].AA3 = (char *)malloc(nNonsy3*sizeof(char));
			G2[dm][k].Sy3 = (char *)malloc(nSynon3*sizeof(char));
			G2[dm][k].AA4 = (char *)malloc(nNonsy4*sizeof(char));
			G2[dm][k].Sy4 = (char *)malloc(nSynon4*sizeof(char));

		}
	}

	Nvrs = (int *)malloc(Ndeme*sizeof(int));

	sTep = ((DinY)*10)/SamInt;
	sampAnces = (double ***)malloc(sTep*sizeof(double **));
	sampEpi = (char ***)malloc(sTep*sizeof(char **));
	for(i = 0; i<sTep; i++){
		sampAnces[i] = (double **)malloc((Ndeme*Nsamp)*sizeof(double *));
		sampEpi[i] = (char **)malloc((Ndeme*Nsamp)*sizeof(char *));
		for(j = 0; j<Ndeme*Nsamp; j++){
			sampAnces[i][j] = (double *)malloc(numEpitope*sizeof(double));
			sampEpi[i][j] = (char *)malloc(numEpitope*sizeof(char));
		}
	}


}

void CleanUp()
{
	int dm, k, loc;
	int sTep = ((DinY)*10)/SamInt;

	for (dm=0; dm<Ndeme; dm++)
	{
		for (k=0; k<Npop; k++)
		{
			free(G1[dm][k].AA);
			free(G1[dm][k].Sy);
			free(G1[dm][k].Rfit);
			free(G1[dm][k].index);

			free(G2[dm][k].AA);
			free(G2[dm][k].Sy);
			free(G2[dm][k].Rfit);
			free(G2[dm][k].index);
		}
		free(G1[dm]);
		free(G2[dm]);
	} 
	free(G1);
	free(G2);

	free (Nvrs);

	for(dm = 0; dm<sTep; dm++){
		for(k=0; k<Ndeme*Nsamp; k++){
			free(sampAnces[dm][k]);
			free(sampEpi[dm][k]);
		}
		free(sampAnces[dm]);
		free(sampEpi[dm]);
	}
	//free(beneAllele);
}

#define IA 16807 
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
        
float ran1( long *idum )
{
        int j;
        long k;
        static long iy = 0;
        static long iv[NTAB];
        float temp;
        
        if (*idum <= 0 || !iy) {
                if (-(*idum) <1 ) *idum =1;
                else *idum = -(*idum);
                for (j=NTAB+7; j>=0; j--) {
                        k= (*idum)/IQ;
                        *idum=IA*(*idum-k*IQ)-IR*k;
                        if (*idum < 0) *idum += IM;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy =iv[0];  
        }
        k=(*idum)/IQ; 
        *idum = IA*(*idum-k*IQ)-IR*k;
        if (*idum < 0) *idum += IM;
        j=iy/NDIV;
        iy=iv[j];
        iv[j] = *idum;
        if ((temp=AM*iy) > RNMX ) return RNMX;
        else return temp;
}

float gammln(float xx)
{
	double x, y, tmp, ser;
	static double cof[6]={76.18009172947146, -86.50532032941677, 24.01409824083091,
						-1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
	int j;
	
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.000000000190015;
	for (j=0; j<=5; j++)
		ser += cof[j]/(++y);
	return -tmp + log(2.5066282746310005*ser/x);
}

#define PI 3.141592654
float poidev( float xm, long *idum )
{
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
		do
		{
			++em;
			t *= ran1(idum);
		}
		while (t > g);
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
		do
		{
			do
			{
				y = tan(PI*ran1(idum));
				em = sq*y+xm;
			}
			while (em < 0.0);
			
			em = floor(em);
			t = 0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		}
		while (ran1(idum) > t);
	}
	return em;
}



void DispVirus(vtype *pvrs)
{
	int loc;

	for (loc=0; loc<nNonsy; loc++)
		fprintf(outfile, "%c", pvrs->AA[loc]);
	fprintf(outfile, " ");
	for (loc=0; loc<nSynon; loc++)
		fprintf(outfile, "%c", pvrs->Sy[loc]);

	for (loc=0; loc<nNonsy; loc++)
		fprintf(outfile, "%c", pvrs->AA2[loc]);
	fprintf(outfile, " ");
	for (loc=0; loc<nSynon; loc++)
		fprintf(outfile, "%c", pvrs->Sy2[loc]);

	for (loc=0; loc<nNonsy3; loc++)
		fprintf(outfile, "%c", pvrs->AA3[loc]);
	fprintf(outfile, " ");
	for (loc=0; loc<nSynon3; loc++)
		fprintf(outfile, "%c", pvrs->Sy3[loc]);

	for (loc=0; loc<nNonsy4; loc++)
		fprintf(outfile, "%c", pvrs->AA4[loc]);
	fprintf(outfile, " ");
	for (loc=0; loc<nSynon4; loc++)
		fprintf(outfile, "%c", pvrs->Sy4[loc]);


	//printf(" %d\n", pvrs->index);
	fprintf(outfile, "\n");
}

void DispPop(vtype **vrs)		//YS//not used
{
	int dm, k, j;

	for (dm=0; dm<Ndeme; dm++)
	{
		printf("deme %d\n", dm+1);
		fprintf(outfile, "deme %d\n", dm+1);
		for (k=0; k<Nvrs[dm]; k++)
			 DispVirus( vrs[dm]+k );
	}
}

float cap(int dm, long tm)		//YS//just returning Kmas
{
	float CC, Day, hSsn;

	Day = (float) (tm%DinY);
	hSsn = DinY*0.33;
/*
	if(dm < Ndeme-2) {
		if(Day < hSsn)
			CC = Kmax*(1.0-Day/hSsn) + 0.01;
		else if(Day > DinY-hSsn)
			CC = Kmax*(Day-DinY+hSsn)/hSsn + 0.01;
		else
			CC = 0.01;
	}else{
		if (Day > DinY/2-hSsn && Day < DinY/2 + hSsn)
			CC = Kmax*(1.0-abs(Day-DinY/2)/hSsn) + 0.01;
		else
			CC = 0.01;
	}
*/
	return Kmax;
}

void Replicate( vtype *parent, vtype *child )	//YS//just copy parent to child
{
	int dm, loc;

	for (loc=0; loc<nNonsy; loc++)
		child->AA[loc] = parent->AA[loc];
	for (loc=0; loc<nSynon; loc++)
		child->Sy[loc] = parent->Sy[loc];

	for (loc=0; loc<nNonsy; loc++)
		child->AA2[loc] = parent->AA2[loc];
	for (loc=0; loc<nSynon; loc++)
		child->Sy2[loc] = parent->Sy2[loc];

	for (loc=0; loc<nNonsy3; loc++)
		child->AA3[loc] = parent->AA3[loc];
	for (loc=0; loc<nSynon3; loc++)
		child->Sy3[loc] = parent->Sy3[loc];

	for (loc=0; loc<nNonsy4; loc++)
		child->AA4[loc] = parent->AA4[loc];
	for (loc=0; loc<nSynon4; loc++)
		child->Sy4[loc] = parent->Sy4[loc];

	for (dm=0; dm<Ndeme; dm++)
		child->Rfit[dm] = parent->Rfit[dm];

	for (loc=0; loc<numEpitope; loc++)
		child->index[loc] = parent->index[loc];

	child->numDel1 = parent->numDel1;
	child->numDel2 = parent->numDel2;
}

void Recombination(vtype *vrs1, vtype *vrs2){		//YS//not used
        int locrc, i;
        unsigned int *temp;

        locrc = (int)(nNonsy*ran1(&seed));
        temp = (unsigned int*)malloc(nNonsy*sizeof(unsigned int));

        for(i=0; i<=locrc; i++){
                temp[i] = vrs1->AA[i];
                vrs1->AA[i] = vrs2->AA[i];
                vrs2->AA[i] = temp[i];
        }
}

int Reassortment(int gen, vtype **vs2, int *nchr){
	int rchr1, rchr2, i, r, numR, dm;
	char temp;
	int tempDel;
	if(doRecomb == 0) { return 0; }
	for (dm=0; dm<Ndeme; dm++){
	numR = poidev(nchr[dm]*0.00001*doRecomb, &seed );

	for (r=0; r< numR; r++){
		rchr1 = (nchr[dm]*ran1(&seed));
		rchr2 = (nchr[dm]*ran1(&seed));

		//first seg
		if(ran1(&seed) < 0.5){
			for(i=0; i<nNonsy; i++){								//YS//exchange rchr1's AA and rchr2's AA
				temp=vs2[dm][rchr1].AA[i];
				vs2[dm][rchr1].AA[i]=vs2[dm][rchr2].AA[i];
				vs2[dm][rchr2].AA[i]=temp;
			}
			tempDel = vs2[dm][rchr1].numDel1;						//YS//exchange rchr1's numDel1 and rchr2's numDel1
			vs2[dm][rchr1].numDel1 = vs2[dm][rchr2].numDel1;
			vs2[dm][rchr2].numDel1 = tempDel;
			for(i=0; i<nSynon; i++){								//YS//exchange rchr1's Sy and rchr2's Sy
				temp=vs2[dm][rchr1].Sy[i];
				vs2[dm][rchr1].Sy[i]=vs2[dm][rchr2].Sy[i];
				vs2[dm][rchr2].Sy[i]=temp;
			}
		}
		//second seg
		if(ran1(&seed) < 0.5){
			for(i=0; i<nNonsy; i++){								//YS//exchange rchr1's AA and rchr2's AA
				temp=vs2[dm][rchr1].AA2[i];
				vs2[dm][rchr1].AA2[i]=vs2[dm][rchr2].AA2[i];
				vs2[dm][rchr2].AA2[i]=temp;
			}
			tempDel = vs2[dm][rchr1].numDel2;						//YS//exchange rchr1's numDel1 and rchr2's numDel1
			vs2[dm][rchr1].numDel2 = vs2[dm][rchr2].numDel2;
			vs2[dm][rchr2].numDel2 = tempDel;
			for(i=0; i<nSynon; i++){
				temp=vs2[dm][rchr1].Sy2[i];							//YS//exchange rchr1's Sy and rchr2's Sy
				vs2[dm][rchr1].Sy2[i]=vs2[dm][rchr2].Sy2[i];
				vs2[dm][rchr2].Sy2[i]=temp;
			}
		}/*
		//third seg
		if(ran1(&seed) < 0.5){
			for(i=0; i<nNonsy3; i++){
				temp=vs2[dm][rchr1].AA3[i];
				vs2[dm][rchr1].AA3[i]=vs2[dm][rchr2].AA3[i];
				vs2[dm][rchr2].AA3[i]=temp;
			}for(i=0; i<nSynon3; i++){
				temp=vs2[dm][rchr1].Sy3[i];
				vs2[dm][rchr1].Sy3[i]=vs2[dm][rchr2].Sy3[i];
				vs2[dm][rchr2].Sy3[i]=temp;
			}
		}
		//fourth seg
		if(ran1(&seed) < 0.5){
			for(i=0; i<nNonsy4; i++){
				temp=vs2[dm][rchr1].AA4[i];
				vs2[dm][rchr1].AA4[i]=vs2[dm][rchr2].AA4[i];
				vs2[dm][rchr2].AA4[i]=temp;
			}for(i=0; i<nSynon4; i++){
				temp=vs2[dm][rchr1].Sy4[i];
				vs2[dm][rchr1].Sy4[i]=vs2[dm][rchr2].Sy4[i];
				vs2[dm][rchr2].Sy4[i]=temp;
			}
		}
		*/
	}
	}
	return numR;
}

void Reproduction(int gen, vtype **vs1, vtype **vs2, int *nchr)
{
        int i, count, dm, j, k;
        float dist, absW, fitn;
        for (dm=0; dm<Ndeme; dm++)
        {
                absW = 2.0/(1+nchr[dm]/cap(dm, gen));
		//printf("cap: %f, absW: %f ", cap(dm, gen), absW); 
                count = 0;
                for (i=0; i<nchr[dm]; i++)
                {
                        fitn = absW*vs1[dm][i].Rfit[dm];
                        j = poidev( fitn, &seed );

                        for (k=0; k<j; k++)
                        {
                                Replicate( vs1[dm]+i, vs2[dm]+count );
                                count++;
                        }
                }

                nchr[dm] = count;
//		printf("nchr[%d]: %d, ", dm, nchr[dm]);
        }
}

void Migration(vtype **vrs, int *nchr)
{
	int j, k, dm, dm2, nmig, loc;
	long mchr;

	for (dm=0; dm<Ndeme; dm++)
	{
		for (dm2=0; dm2<Ndeme; dm2++)
		{
			if (dm != dm2)
			{
				nmig = (int) poidev( r_mig[dm][dm2]*nchr[dm], &seed );
				while( nmig>nchr[dm]){
					nmig--;
				}

				for (j=0; j<nmig; j++)
				{
					mchr = (int) (nchr[dm]*ran1(&seed));

					for (loc=0; loc<nNonsy; loc++)
						vrs[dm2][nchr[dm2]].AA[loc] = vrs[dm][mchr].AA[loc];

					for (loc=0; loc<nSynon; loc++)
						vrs[dm2][nchr[dm2]].Sy[loc] = vrs[dm][mchr].Sy[loc];

					for (k=0; k<Ndeme; k++)
						vrs[dm2][nchr[dm2]].Rfit[k] = vrs[dm][mchr].Rfit[k];

					for(loc=0; loc<numEpitope; loc++)
						vrs[dm2][nchr[dm2]].index[loc] = vrs[dm][mchr].index[loc];

					nchr[dm2]++;

					for (loc=0; loc<nNonsy; loc++)
						vrs[dm][mchr].AA[loc] = vrs[dm][nchr[dm]-1].AA[loc];

					for (loc=0; loc<nSynon; loc++)
						vrs[dm][mchr].Sy[loc] = vrs[dm][nchr[dm]-1].Sy[loc];

					for (k=0; k<Ndeme; k++)
						vrs[dm][mchr].Rfit[k] = vrs[dm][nchr[dm]-1].Rfit[k];

					for(loc=0; loc<numEpitope; loc++)
						vrs[dm][mchr].index[loc] = vrs[dm][nchr[dm]-1].index[loc];

					nchr[dm]--;
				}
			}
		}
	}

}	

void Sampling(int tm, vtype **vrs, int *nchr, double ***sA, char ***sE, int *extinct)
{
	int dm, j, k, nsam, vrsp, stime, i;
	int schr, Gnchr, sdm;

	vrsp = 0;
	j = 0;
	while ( vrsp < 1 && j < Ndeme )
	{
		if ( nchr[j++] > 0 )
			vrsp = 1;
	}

	if (vrsp == 0)
	{
		*extinct = 1;
	}
	else
	{
		if (tm%SamInt==0)
		{
			stime = (tm%(2*(DinY*10)))/SamInt;
			Gnchr = nchr[0];

			for (dm=1; dm<Ndeme; dm++)
			{
				if(nchr[dm] > Gnchr){
					Gnchr = nchr[dm];
				}
			}

			for(j=0; j<Ndeme*Nsamp; j++)
			{
				sdm = (int) (Ndeme*ran1(&seed));
				while(1.0*nchr[sdm]/Gnchr < ran1(&seed)){
					sdm = (int)(Ndeme*ran1(&seed));
				}

				schr = (int) (nchr[sdm]*ran1(&seed));
/*				for(i=0; i<numEpitope; i++){
					sA[stime][j][i] = vrs[sdm][schr].index[i];
					sE[stime][j][i] = vrs[sdm][schr].AA[i];
					//printf("gen: %d, sA[%d][%d][%d] : %lf, sE : %c ", tm+Tburn, stime, j, i, sA[stime][j][i], sE[stime][j][i]);
				}*/
				if (dispv > 0 )
				{
					//printf("g%d d%d\n", (stime+1)*SamInt, dm+1);
					fprintf(outfile, "g%d d%d %lf %lf %lf ", tm%(2*DinY*10), sdm+1, vrs[sdm][schr].index[0], vrs[sdm][schr].index[1], vrs[sdm][schr].index[2]);
					DispVirus( &vrs[sdm][schr] );
				}
			}//end of for j

		}//end of if tm
	}
}

void Mutation(vtype **vrs, int *nchr, int thisgen)
{
	int dm, nmut, i, j, mchr, mloc, mUnit, episite;
	unsigned int seq, mutint;

	for (dm=0; dm<Ndeme; dm++)
	{
		nmut = (int) poidev(r_mut*nNonsy*nchr[dm], &seed);

		for (j=0; j<nmut; j++)
		{
			mloc = (int) (nNonsy*ran1(&seed));
			episite = 0;
			if(mloc < numEpitope){
				if( (int)(reduced_u*ran1(&seed)) != 0 ){
						episite = 1;
						continue;
				}
			}

			mchr = (int) (nchr[dm]*ran1(&seed));

			if(vrs[dm][mchr].AA[mloc] == '0') {
				if( mloc >= numEpitope || ( mloc < numEpitope && vrs[dm][mchr].AA[mloc] != beneAllele[mloc])){
					vrs[dm][mchr].AA[mloc] = '1';
				}
				if (mloc >= numEpitope) {
					vrs[dm][mchr].numDel1 += 1;
				}
			}else{
				if( mloc >= numEpitope || ( mloc < numEpitope && vrs[dm][mchr].AA[mloc] != beneAllele[mloc])){
					vrs[dm][mchr].AA[mloc] = '0';
				}
				if (mloc >= numEpitope) {
					vrs[dm][mchr].numDel1 -= 1;
				}
			}

		}

		nmut = (int) poidev(r_mut*nSynon*nchr[dm], &seed);

		for (j=0; j<nmut; j++)
		{
			mloc = (int) (nSynon*ran1(&seed));
			mchr = (int) (nchr[dm]*ran1(&seed));

			if(vrs[dm][mchr].Sy[mloc] == '0') {
				vrs[dm][mchr].Sy[mloc] = '1';
			}else{
				vrs[dm][mchr].Sy[mloc] = '0';
			}
		}

		//second segment
		nmut = (int)poidev(r_mut*nNonsy*nchr[dm], &seed);
		for(j=0; j<nmut; j++){
			mloc = (int)(nNonsy*ran1(&seed));
			mchr = (int)(nchr[dm]*ran1(&seed));
			if(vrs[dm][mchr].AA2[mloc] == '0'){
				vrs[dm][mchr].AA2[mloc]='1';
				vrs[dm][mchr].numDel2 += 1;
			}
			else{
				vrs[dm][mchr].AA2[mloc]='0';
				vrs[dm][mchr].numDel2 -= 1;
			}
		}
		nmut = (int)poidev(r_mut*nSynon*nchr[dm], &seed);
		for(j=0; j<nmut; j++){
			mloc = (int)(nSynon*ran1(&seed));
			mchr = (int)(nchr[dm]*ran1(&seed));
			if(vrs[dm][mchr].Sy2[mloc] == '0'){vrs[dm][mchr].Sy2[mloc]='1';}
			else{vrs[dm][mchr].Sy2[mloc]='0';}
		}
		/*
		//third segment
		nmut = (int)poidev(r_mut*nNonsy3*nchr[dm], &seed);
		for(j=0; j<nmut; j++){
			mloc = (int)(nNonsy3*ran1(&seed));
			mchr = (int)(nchr[dm]*ran1(&seed));
			if(vrs[dm][mchr].AA3[mloc] == '0'){
				vrs[dm][mchr].AA3[mloc]='1';
				vrs[dm][mchr].numDel += 1;
			}
			else{
				vrs[dm][mchr].AA3[mloc]='0';
				vrs[dm][mchr].numDel -= 1;
			}
		}
		nmut = (int)poidev(r_mut*nSynon3*nchr[dm], &seed);
		for(j=0; j<nmut; j++){
			mloc = (int)(nSynon3*ran1(&seed));
			mchr = (int)(nchr[dm]*ran1(&seed));
			if(vrs[dm][mchr].Sy3[mloc] == '0'){vrs[dm][mchr].Sy3[mloc]='1';}
			else{vrs[dm][mchr].Sy3[mloc]='0';}
		}
		//fourth segment
		nmut = (int)poidev(r_mut*nNonsy4*nchr[dm], &seed);
		for(j=0; j<nmut; j++){
			mloc = (int)(nNonsy4*ran1(&seed));
			mchr = (int)(nchr[dm]*ran1(&seed));
			if(vrs[dm][mchr].AA4[mloc] == '0'){
				vrs[dm][mchr].AA4[mloc]='1';
				vrs[dm][mchr].numDel += 1;
			}
			else{
				vrs[dm][mchr].AA4[mloc]='0';
				vrs[dm][mchr].numDel -= 1;
			}
		}
		nmut = (int)poidev(r_mut*nSynon4*nchr[dm], &seed);
		for(j=0; j<nmut; j++){
			mloc = (int)(nSynon4*ran1(&seed));
			mchr = (int)(nchr[dm]*ran1(&seed));
			if(vrs[dm][mchr].Sy4[mloc] == '0'){vrs[dm][mchr].Sy4[mloc]='1';}
			else{vrs[dm][mchr].Sy4[mloc]='0';}
		}
	*/
	}
}

void Initialize()
{
	int dm, k, loc, j;

	for (dm=0; dm<Ndeme; dm++)
	{
		Nvrs[dm] = (int) Kmax;

		for (k=0; k<Nvrs[dm]; k++)
		{
			for (loc=0; loc<nNonsy; loc++)
				G1[dm][k].AA[loc] = '0';
			for (loc=0; loc<nSynon; loc++)
				G1[dm][k].Sy[loc] = '0';
			for (loc=0; loc<nNonsy; loc++)
				G1[dm][k].AA2[loc] = '0';
			for (loc=0; loc<nSynon; loc++)
				G1[dm][k].Sy2[loc] = '0';
			for (loc=0; loc<nNonsy3; loc++)
				G1[dm][k].AA3[loc] = '0';
			for (loc=0; loc<nSynon3; loc++)
				G1[dm][k].Sy3[loc] = '0';
			for (loc=0; loc<nNonsy4; loc++)
				G1[dm][k].AA4[loc] = '0';
			for (loc=0; loc<nSynon4; loc++)
				G1[dm][k].Sy4[loc] = '0';

			for (j=0; j<Ndeme; j++)
				G1[dm][k].Rfit[j] = 1.0;
			for (j=0; j<numEpitope; j++)
				G1[dm][k].index[j] = 0.0;
				
			G1[dm][k].numDel1 = 0;
			G1[dm][k].numDel2 = 0;
		}
	}

}



void setBeneAllele() {
	int i;
	for(i=0; i<numEpitope; i++) {
		beneAllele[i] = '1';
	}
}

void UpdateFitness(vtype **vrs, int *nchr){

	int i, j, d, dm, fUnit, flocation, numAllVrs;
	unsigned int seq, fitint;
	float sumFitness, meanFitness;
	mostfit = 1.0;
	//initialize the fitness
	for(dm=0; dm<Ndeme; dm++){
		for(j=0; j<nchr[dm]; j++){
			for(d=0; d<Ndeme; d++){
				vrs[dm][j].Rfit[d] = 1.0;
			}
		}
	}

//		printf("floc: %d\n", floc[i]);
	sumFitness = 0;
	meanFitness = 0;
	numAllVrs = 0;
	for (dm=0; dm<Ndeme; dm++)
	{
		numAllVrs += nchr[dm];
		for (j=0; j<nchr[dm]; j++)
		{
			for(i=0; i<numEpitope; i++) {
				if (vrs[dm][j].AA[i] == beneAllele[i]){
					for(d=0; d<Ndeme; d++){
						vrs[dm][j].Rfit[d] *= (1+sel_co);
					}
				}
			}

			for(d=0; d<Ndeme; d++){
					vrs[dm][j].Rfit[d] *= pow( (1-sel_d), (vrs[dm][j].numDel1 + vrs[dm][j].numDel2));
			}

			sumFitness += vrs[dm][j].Rfit[dm];
		}//end for j individual
	}

	meanFitness = sumFitness/(1.0*numAllVrs);

	for(dm=0; dm<Ndeme; dm++){
		for(j=0; j<nchr[dm]; j++){
			for(d=0; d<Ndeme; d++){
				vrs[dm][j].Rfit[d] = vrs[dm][j].Rfit[d]/meanFitness;
				if (vrs[dm][j].Rfit[d] > mostfit) {
					mostfit = vrs[dm][j].Rfit[d];
				}
			}
		}
	}

}

int softsweep(double ***sA, char ***sE, int site, int gen){
	int Laststime, earlyC, Cstime;
	int j, g, l, countA, countC, numlin, found;
	int linfreq[100];
	double linhis[100];
	int majorfreq;

	//get Laststime
	Laststime = ((gen-Tburn)%(2*DinY*10))/SamInt;
	if (Laststime <= 0 || Laststime >= (DinY*10)/SamInt){
		printf("out of observation at gen: %d\n", gen);
		return 0;
	}

	//get Cstime
	earlyC = gen;
	for(j=0; j<Ndeme*Nsamp; j++){
		if(sE[Laststime][j][site] == '1' && (sA[Laststime][j][site]) < earlyC){
			earlyC = (int)(sA[Laststime][j][site]);
		}
	}
	if( (earlyC-Tburn)<0 || ((earlyC-Tburn)/(2*DinY*10) != (gen-Tburn)/(2*DinY*10)) ){
		printf("out of observation at gen: %d\n", gen);
		return 0;
	}
	Cstime = ((earlyC-Tburn)%(2*DinY*10))/SamInt;
	printf("checking ss at site: %d, gen: %d, try: %d, Laststime: %d, Cstime: %d\n", site, gen, 1+(gen-Tburn)/(2*DinY*10), Laststime, Cstime);

	numlin = 0;
	for(g = Cstime; g <= Laststime; g++){
		countA = 0;
		countC = 0;
		for(l=0; l<numlin; l++){
			linfreq[l] = 0;
		}
		for(j=0; j<Ndeme*Nsamp; j++){
			if (sE[g][j][site] == '0') {
				countA += 1;
			}else {
				countC += 1;
				if(numlin==0){
					linhis[0] = sA[g][j][site];
					linfreq[0] = 1;
					numlin = 1;
				}else{
					found = 0;
					for(l=0; l<numlin; l++){
						if(sA[g][j][site] == linhis[l]){
							linfreq[l] += 1;
							found = 1;
							break;
						}
					}
					if(found==0){
						linhis[numlin] = sA[g][j][site];
						linfreq[numlin] = 1;
						numlin += 1;
					}
									}
			}
		}// end of for j

		printf("ptime: %d, stime: %d, A: %d, C: %d | ", g*6, g, countA, countC); 
		majorfreq = 0;
		for(l=0; l<numlin; l++){
			if(linfreq[l] == 0){
				continue;
			}
			printf("linhis[%d]: %lf, freq: %d | ", l, linhis[l], linfreq[l]);
			if(linfreq[l] > majorfreq) {
				majorfreq = linfreq[l];
			}
		}
		printf("\n");
/*		printf("ptime: %d, stime: %d, countA: %d, countC: %d\n", g*6, g, countA, countC); 
		if(1.0*countC/(countC+countA) >= 0.5 && 1.0*majorfreq/countC < 0.8){
			printf("softsweep! site: %d, C freq = %lf, major lineage freq = %lf\n", site, 1.0*countC/(countA+countC), 1.0*majorfreq/countC);
			break;
		}
*/
	}//end of for g
	return Laststime;

}

int ifFixed_resetEpisite(vtype **vrs, int *nchr, double ***sA, char ***sE, int gen){
	int i, j, dm, e, nchrAll; //for iteration
	int count, numUpdated;
	int *fixed = (int*)malloc(nNonsy*sizeof(int));
	unsigned int seq, fitint;

	numUpdated = 0;

	for(i=0; i<numEpitope; i++) {
		fixed[i] = 0;
		count = 0;
		nchrAll = 0;

		for (dm=0; dm<Ndeme; dm++)
		{
			nchrAll += nchr[dm];
			for (j=0; j<nchr[dm]; j++)
			{
				if(vrs[dm][j].AA[i] == beneAllele[i]) {
					count += 1;
				}
			}
		}
		if(count < nchrAll) { //not fixed
			continue;
		}
		else {
			fixed[i] = 1;
			numUpdated+=1;
			if (beneAllele[i] == '1') { beneAllele[i] = '0';}
			else { beneAllele[i] = '1';} 
			//softsweep(sA, sE, i, gen);

		}//end of else
	}//end of for numEpi
	if(numUpdated>=1) {
		printf("gen: %d , numUpdated: %d |  ", gen, numUpdated);
		for(e=0; e<numEpitope; e++) {
			if (fixed[e] == 1){
				printf("%d ", e);
			} 
		}
		printf("\n");
	}
	return numUpdated;
}

main(int argc, char*argv[])
{
	int i, j, e, extc, gen, nseq, Tnew, numUpdated, avgSubinEpi;
	int vrsp, dm, extincted, realNtry; // for checking extinction
	float fstN, fstS, cdist, pi;
	long thisSeed;
	const char* inpname;
	const char* outname;
	vtype **temp, **vrs1, **vrs2;
	printf("main\n");
	if(argc==4){
		thisSeed = atoi(argv[3]);
		inpname = argv[1];
		outname = argv[2];
		printf("%s\n", argv[1]);
	}
	else if (argc==1) {
		thisSeed = 0;
		inpname = "fsim.inp";
		outname = "fsim.out";
	}
	Parameters(inpname, outname, thisSeed);
	GetArrays();
	avgSubinEpi = 0;
	realNtry = 0;

	for (i=0; i<Ntry; i++)
	{
		extincted = 0;
		printf("\nTry%d ", i+1);
		//fprintf(outfile, "\n*\nTry%d\n", i+1);
		vrs1 = G1;
		vrs2 = G2;
		Initialize();
		setBeneAllele();
		gen = 0;
		numUpdated = 0;

		printf("\n");

		extc = 0;
		while( gen < Tmax && !extc )
		{

			vrsp = 0;
			dm = 0;
			while ( vrsp < 1 && dm < Ndeme )
			{
				if ( Nvrs[dm++] > 0 )
					vrsp = 1;
			}
			if (vrsp == 0) {
				printf("extinct!!\n");
				extincted = 1;
				break;
			}

//			Migration(vrs1, Nvrs);
			Mutation(vrs1, Nvrs, gen);
			Reproduction(gen, vrs1, vrs2, Nvrs);
			Reassortment(gen, vrs2, Nvrs);
			gen++;
			temp = vrs1;
			vrs1 = vrs2;
			vrs2 = temp;
			if (gen >= Tburn){
				if((gen-Tburn)%(2*DinY*10) == 0){
					fprintf(outfile, "\n*\nTry%d\n", 1+(gen-Tburn)/(2*DinY*10));
				}
				if((gen-Tburn)%(2*DinY*10) < DinY*10){
					Sampling(gen-Tburn, vrs1, Nvrs, sampAnces, sampEpi, &extc);
				}
			}

			numUpdated += ifFixed_resetEpisite(vrs1, Nvrs, sampAnces, sampEpi, gen);
			UpdateFitness(vrs1, Nvrs);

			if ( gen%DinY == dWin/2 && gen >= Tburn)
			{
				Tnew = gen-Tburn;

				printf("%d %d ", Tnew, extc);
				for (j=0; j<Ndeme; j++)
					printf(" %5.0f %5d", cap(j, gen), Nvrs[j]);
				printf("\n");

				//CalcHetero(2, &pi, &nseq, -1, Tnew-dWin, Tnew);
				//printf("Pi= %f, Nseq= %d, ", pi, nseq);
			
				//fstN = CalculateFst(1, Tnew-dWin, Tnew);
				//fstS = CalculateFst(0, Tnew-dWin, Tnew);
				//printf("FstN= %f, FstS= %f, ratio=%f \n", fstN, fstS, fstN/fstS);

				//cdist = CalculateDst(2, gen, vrs1, Nvrs);
				//printf("cdist = %f mut*gen = %f\n", cdist, r_mut*gen);
			}

			//if ( gen >= Tburn && gen%DinY == dWin/2 )
				//GetHad( gen-Tburn-dWin, gen-Tburn);
		}
		if(extincted==0){
			printf("******* number of susbstitution in epitopes: %d\n\n", numUpdated); 
			avgSubinEpi += numUpdated;
			realNtry += 1;
		}
	}//end of for try
	avgSubinEpi /= realNtry;
	printf("\n realNtry: %d, avgSubinEpi: %d\n", realNtry, avgSubinEpi);
	//CalcMeans();

	CleanUp();
	return 0;
}
