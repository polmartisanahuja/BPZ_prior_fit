#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
//#include "/Users/pmarti/Dropbox/Pol/Tesi/photoz/2slaq_DR7_LRG_webcat/BPZ/PRIOR/prior_fit.h"
//include "../PhD_pmarti/Routine/Parameters/pau_prior_fit.h"
#include "../Parameters/pau_prior_fit.h"

/*This routine "prior_fit_v1" works together with "nrutil.c". It fits the prior function described in Benitez 2000 using the values: 
magnitude, redshift and spectral type of a collection of galaxies. To do this, a likelihood is maximized through
a method called Direction set (POWELL's) Meethod also used in the prior calibration did in Benitez 2000. POWELL's meethod is obtained 
and copied from the Numerical Recipies. The routine also generates a plot of the fitted curve against the true redshift and type distribution */

#define ITMAX 200 //Maximum allowed iterations.
#define TOL 2.0e-4 //Tolerance passed to brent.

FILE *file_init(char file_name[50]);
int proportions();
int reductor(int N);
float function(float *lam);
void powell(float p[], float **xi, int n, float ftol, int *iter, float *fret, float (*func)(float []));
void printprior(float *lam);
void plotprior(float *lam);

//Global variables
float *m_smallcat, *z_smallcat, *t_smallcat;
int N_smallcat;
double f[3];
float f_min = 10000000000;
float ftol = 0.001; //(Default ftol=0.000001) It can not be set as a constant (#define)
int n_para = 11; //Number of parameters of the prior

main()
{
	FILE *priorpara, *mock;
	char word[10], cmd[50];
	int N, i, j;
	int *iter;
	float *lam, **xi, *fret;
	
	printf("\n***********************************************\n");
	printf("*                                             *\n");
	printf("*       PRIOR CALIBRATION by P. Marti         *\n");
	printf("*                 Version 1                   *\n");
	printf("*                                             *");
	printf("\n***********************************************\n");
	
	//Reserving memory for lam[]....................
	lam = (float*) calloc(n_para+1,sizeof(float));
		
	if(display != 0){
	
		//Read Prior parameters file...........................
		priorpara = file_init(prior_para_file);
		for(i = 1; i <= 2; i++) fscanf(priorpara,"%s %lf\n", word, &f[i]);
		for(i = 1; i <= n_para; i++) fscanf(priorpara,"%s %f\n", word, &lam[i]);
		fclose(priorpara);
	}
	else{
		
		//Proportions..........................................
		N = proportions();
		
		if(N <= n_gal_fit){
			printf("There isn't enought galaxies to do the fit.\n");
			exit(1);}
		
		//Reductor and m cut...................................
		N_smallcat = reductor(N);

		//Reserving memory and seeting the m,z,t small_cat values..........................
		i = 0;
		mock = file_init(small_cat_file);
		m_smallcat = (float*) calloc(N_smallcat,sizeof(float));
		z_smallcat = (float*) calloc(N_smallcat,sizeof(float));
		t_smallcat = (float*) calloc(N_smallcat,sizeof(float));
		while(fscanf(mock, "%f %f %f\n", &m_smallcat[i], &z_smallcat[i], &t_smallcat[i]) != -1) i++;
		fclose(mock);

		//Reserving memory and setting xi, lam and frt, iter variables...............
		xi = (float**) calloc(n_para+1,sizeof(float));
		for(i = 0; i <= n_para; i++) xi[i] = (float*) calloc(n_para+1,sizeof(float));

		for(i = 1; i <= n_para; i++){
			for(j = 1; j <= n_para; j++){
				if(j == i) xi[j][i] = 0.1;
				else xi[j][i] = 0;}}
		
		//Seeting lam[].........................................
		for(i = 1; i <= n_para; i++) lam[i] = 0.1;
		
		//Reserving memory for fret and iter....................
		fret = (float*) calloc(1,sizeof(float));
		iter = (int*) calloc(1,sizeof(int));
		
		//Minimizing the function...............................
		printf("\n---MINIMIZING MODULE----------------------------------------\n");
		printf("\tComputing...\n");
		powell(lam,xi,n_para,ftol,iter,fret, function);
		
		printf("\n\tfret=%f\n", fret[0]);
		//printf("\niter=%f\n", iter[0]);
			
		//Make absolute values of the parameters................
		for(i = 1;i <= n_para;i++) lam[i] = fabs(lam[i]);
		
		priorpara=fopen(prior_para_file,"w");
		fprintf(priorpara,"#Fitted function parameters\n");
		for(i = 1;i <= 2;i++) fprintf(priorpara,"f_%d\t%3.3lf\n", i, f[i]);
		for(i = 1;i <= 2;i++) fprintf(priorpara,"k_%d\t%3.3f\n", i, lam[i]);
		for(i = 1;i <= 3;i++) fprintf(priorpara,"a_%d\t%3.3f\n", i, lam[i+2]);
		for(i = 1;i <= 3;i++) fprintf(priorpara,"z0_%d\t%3.3f\n", i, lam[i+5]);
		for(i = 1;i <= 3;i++) fprintf(priorpara,"k_m%d\t%3.3f\n", i, lam[i+8]);
		fclose(priorpara);
	}
		
	//Write Prior parameters in the terminal and in a file...........................
	printf("\n\tFitted function parameters:\n");
	for(i = 1;i <= 2;i++) printf("\tf[%d]\tf_%d\t%3.3lf\n", i, i, f[i]);
	for(i = 1;i <= 2;i++) printf("\tlam[%d]\tk_%d\t%3.3f\n", i, i, lam[i]);
	for(i = 1;i <= 3;i++) printf("\tlam[%d]\ta_%d\t%3.3f\n", i+2, i, lam[i+2]);
	for(i = 1;i <= 3;i++) printf("\tlam[%d]\tz0_%d\t%3.3f\n", i+5, i, lam[i+5]);
	for(i = 1;i <= 3;i++) printf("\tlam[%d]\tk_m%d\t%3.3f\n", i+8, i, lam[i+8]);
	
	//Printing the prior function as .py script for BPZ..............................
	printprior(lam);
	
	//Plotting the prior.............................................................
	plotprior(lam);
	
	//sprintf(cmd, "rm %sp_z_t_*.ps", plot_file);
	//system(cmd);

return 0;
}



//FUNCTIONS.....................................................................................................................

FILE *file_init(char file_name[50])
{
	/*Returns a pointer to the file "in" which has been initialized to the first 
	row where data is. This means that all headers of the file that start with # 
	are avoided*/
	
	FILE *in;
	int N = 0, i = 0;
	char caracter1, caracter2;
	
	//Count number of lines begining with # (headers)
	if ((in=fopen(file_name, "r")) == NULL) printf ("Error opening file \"%s\".\n",file_name);

	fscanf(in, "%c", &caracter1);
	while(caracter1=='#') {
		N++;
		fscanf(in, "%c", &caracter2);
		while(caracter2!='\n') fscanf(in, "%c", &caracter2);
		fscanf(in, "%c", &caracter1);}
	fclose(in);
	
	//Read rows until data will be found
	if ((in=fopen(file_name, "r")) == NULL) printf ("Error opening file %s.\n",file_name);

	while(i != N) {
		fscanf(in, "%c", &caracter1);
		if(caracter1 == '\n') i++;}
		
	return in;
}

int proportions()
{
	/*This function returns the total number of objects in the catalog and also compute
	the fraction f(m=m_0) of early type galaxies as well as the Spirals, which will become
	parameters of the prior function. This function also prints on the terminal useful
	information as the magnitude and type range of the catalog.*/

	FILE *mock;
	double m, z, t, m_min, m_max, t_min, t_max, Delta = 0.15;
	int error, n = 0, N = 0;
	
	//Find f[1] and f[2].................
	f[1]=0;
	f[2]=0;
	mock = file_init(cat_file);
	error = fscanf(mock, "%lf %lf %lf\n", &m, &z, &t);
	m_min = m;
	m_max = m;
	t_min = t;
	t_max = t;
	while(error != -1){
		if( (m < (m_0+Delta)) && (m > (m_0-Delta)) ){
			if(t <= t1) f[1]++;
			else if(t <= t2) f[2]++;
			n++;}
		N++;
		if(m_min > m) m_min = m;
		if(m_max < m) m_max = m;
		if(t_min > t) t_min = t;
		if(t_max < t) t_max = t;
		//printf("N=%d\n", N);
		error = fscanf(mock, "%lf %lf %lf\n", &m, &z, &t);}
	fclose(mock);
	
	f[1] = f[1]/n;
	f[2] = f[2]/n;
	
	//Printing results......................
	printf("\n---PROPORTIONS MODULE----------------------------------------\n");
	printf("\tTotal N gal. in the cat.:\t\t%d\n\tN gal. (at m=%2.2f):\t\t\t%d {m|%2.2lf < m < %2.2lf}\n", N, m_0, n, m_0-Delta, m_0+Delta);
	printf("\tTotal mag. rang:\t\t\t[%2.2lf,%2.2lf]\n", m_min, m_max);
	printf("\tTotal type. rang:\t\t\t[%2.2lf,%2.2lf]\n\n", t_min, t_max);
	printf("\tEarly type gal. prop. (at m=%2.2f):\t%2.2lf%% {t|%2.2lf < t <= %2.2f}\n", m_0, f[1]*100, t_min, t1);
	printf("\tSpiral gal. prop. (at m=%2.2f):\t\t%2.2lf%% {t|%2.2lf < t <= %2.2f}\n\n", m_0, f[2]*100, t1, t2);
	return N;
}

int reductor(int N)
{
	/*This function randomly selects "n_gal_fit" galaxies that will be
	used for the prior fitting. This galaxies are in the magnitude range
	[m_0, m_cut]. It returns this value and writes the new
	reduced catalog in a file called "mzt_small.txt"*/

	FILE *mock, *out;
	int i, n = 0;
	float *m, *z, *t;
	
	m = (float*) calloc(N, sizeof(float));
	z = (float*) calloc(N, sizeof(float));
	t = (float*) calloc(N, sizeof(float));

	mock = file_init(cat_file);
	for(i=0; i<N; i++) fscanf(mock, "%f %f %f\n", &m[i], &z[i], &t[i]);
	fclose(mock);

	out = fopen(small_cat_file,"w");
	fprintf(out, "#%f<m_prior<%f\n#m_prior z_true t_true\n", m_0, m_cut);
	while(n != n_gal_fit){
		i = rand() % N;
		if((m[i]>m_0) && (m[i]<m_cut)){ //Magnitude Cuts
			n++;
			fprintf(out, "%f %f %f\n", m[i], z[i], t[i]);
		}
	}
	fclose(out);

	printf("\n---REDUCTOR MODULE------------------------------------------\n");
	printf("\tTotal N gal. in the cat.:\t\t%d\n\tRandomly selected gal. for prior:\t%d\n\n", N,n);
	
	return n;
}

void printprior(float *lam)
{
	/*This function prints a file.py with the prior function that BPZ
	will use.*/
	
	FILE *prior;
	char cmd[1000];
	
	printf("\n---PRINTING PRIOR.py MODULE--------------------------------\n");
	
	prior = fopen(prior_py_file,"w");
	
	//bpz 1.99.3
	fprintf(prior,"a = %3.3f, %3.3f, %3.3f\n", lam[3], lam[4], lam[5]);
	fprintf(prior,"zo = %3.3f, %3.3f, %3.3f\n", lam[6], lam[7], lam[8]);
	fprintf(prior,"km = %3.3f, %3.3f, %3.3f\n", lam[9], lam[10], lam[11]);
	fprintf(prior,"k_t = %3.3f, %3.3f\n", lam[1], lam[2]);
	fprintf(prior,"fo_t = %2.2f, %2.2f\n\n", f[1], f[2]/2);

	//bpz 1.98b
	//fprintf(prior,"from bpz_tools import *\ndef function(z,m,nt):\n\n\tmmax=28.\n\n\tif nt<>6:\n\t\tprint \"Wrong number of template spectra!\"\n\t\tsys.exit()\n\n\tglobal zt_at_a\n\tnz=len(z)\n");
	//fprintf(prior,"\tmomin_hdf=%2.2f\n", m_0);
	//fprintf(prior,"\tif m>32.: m=32.\n\tif m<%2.2f: m=%2.2f\n", m_0, m_0);
	fprintf(prior,"\ta=array((%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f))\n", lam[3], lam[4], lam[4], lam[5], lam[5], lam[5]);
	fprintf(prior,"\tzo=array((%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f))\n", lam[6], lam[7], lam[7], lam[8], lam[8], lam[8]);
	fprintf(prior,"\tkm=array((%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f))\n", lam[9], lam[10], lam[10], lam[11], lam[11], lam[11]);
	fprintf(prior,"\tfo_t=array((%2.2f,%2.2f,%2.2f))\n", f[1], f[2]/2, f[2]/2);
	fprintf(prior,"\tk_t=array((%3.3f,%3.3f,%3.3f))\n", lam[1], lam[2], lam[2]);
	//fprintf(prior,"\tdm=m-momin_hdf\n\tzmt=clip(zo+km*dm,0.01,15.)\n\tzmt_at_a=zmt**(a)\n");
	//fprintf(prior,"\n\ttry:\n\t\tzt_at_a.shape\n\texcept NameError:\n\t\tzt_at_a=power.outer(z,a)\n");
	//fprintf(prior,"\n\tf_t=zeros((len(a),),Float)\n\tf_t[:3]=fo_t*exp(-k_t*dm)\n\tf_t[3:]=(1.-add.reduce(f_t[:3]))/3.\n");
	//fprintf(prior,"\n\n\tp_i=zt_at_a[:nz,:6]*exp(-clip(zt_at_a[:nz,:6]/zmt_at_a[:6],0.,700.))\n");
	//fprintf(prior,"\n\tnorm=add.reduce(p_i[:nz,:6],0)\n\tp_i[:nz,:6]=p_i[:nz,:6]/norm[:6]*f_t[:6]\n\treturn p_i\n");
	
	fclose(prior);
	sprintf(cmd,"cat %s", prior_py_file);
	system(cmd);
}


float function(float *lam_fit)
{
	/*This is the function that is going to be minimized. In this case is the bayesian function
	described in Benitez 2000 that works as a prior in the BPZ code. "lam_fit" are the 11 free 
	paremeters of the function*/
	
	float f_out=0;
	float pt, pz;
	float zi, zmt, dz, sum;
	float dm;
	int error, t, i, j;
	
	//computing p(z,t|m)=p(t|m)*p(z|t,m)...............................
	for(i = 0; i < N_smallcat; i++){
		dm = m_smallcat[i] - m_0;
		
		if(dm > 0){ //In fact this is not necessary. This cut has been already done in reduct function.
			
			//Setting t value..........................................
			if(t_smallcat[i] <= t1) t = 1;
			else{
				if(t_smallcat[i] <= t2) t = 2;
				else t=3;}
			
			//Forcing parameters to be possitive......................
			for(j = 1;j <= 11; j++) lam_fit[j] = fabs(lam_fit[j]);
			
			//Computing p(t|m).........................................
			if(t<3) pt=f[t]*exp(-lam_fit[t]*dm);
			else pt=1-(f[1]*exp(-lam_fit[1]*dm))-(f[2]*exp(-lam_fit[2]*dm));
	
			//Computing p(z|t,m)......................................
			zmt=lam_fit[t+5]+(lam_fit[t+8]*dm);
			
			//Normalization
			sum=0, zi = 0;
			dz=0.01;
			while(zi <= 20.0){ //20.0 is an enough large redshift 
				sum += (dz*pow(zi,lam_fit[t+2])*exp(-pow(zi/zmt, lam_fit[t+2])));
				zi += dz;}
			pz=(1/sum)*(pow(z_smallcat[i],lam_fit[t+2])*exp(-pow(z_smallcat[i]/zmt,lam_fit[t+2]))); 
			
			f_out+=-log(pt)-log(pz);
			//printf("pt=%f sum=%f pz=%f f_out=%f\n\n",pt,sum,pz, f_out);
		}
	}
	
	if(f_out < f_min) f_min = f_out;
	printf("\tf_min=%3.3f\r", f_min);
	fflush(stdout);
	
	/*//Asertion statment.................................................
	if(fabs(f_out - 2050.284180) < 0.00001) 
	{
		printf("Correct!\n");
		exit(0);
	}*/
	
	return f_out;
}

void plotprior(float *lam)
{
	//This function prints a lis of plots with the prior function obtained from current data,
	//fitted parameters and default Txitxo parameters. Different plots correspond to different magnitudes.
	//Each plot contains 4 panels, one for the type distribution and the other 3 for the redshift distribution
	//for each type.
	
	FILE *mock;
	FILE *p_t_m_file;
	FILE *p_z_t_m_file;
	FILE *p_t_m_func_file;
	FILE *p_t_m_func_tx_file;
	FILE *p_z_t_m_func_file;
	FILE *p_z_t_m_func_tx_file;
	FILE *plot;
		
	float m, z, t, m_max, zmax, t_max, m_min, zmin, t_min, error, delta_z, delta_m, delta_t, norm, histo[bins_m][bins_z][bins_t], p[3], p_tm[2], p_ztm, m_prior;
	float lam_tx[12], f_tx[3], m_0_tx, dm, dm_tx, zi, sum, sum_tx, dz, zmt, zmt_tx;
	float origin[3][2];
	
	int i, j, k;
	int l, n_fora, *bins, N, t_prior;
	
	char nom_fich_temp[1000];
	char nom_fich_tm[1000], nom_fich_tm_func[1000], nom_fich_tm_func_tx[1000];
	char nom_fich_ztm[1000], nom_fich_ztm_func[1000], nom_fich_ztm_func_tx[1000];

	//Set txitxo prior values..................................
	m_0_tx = 20.0;
	f_tx[1] = 0.35, f_tx[2] = 0.5; //f_1, f_2
	lam_tx[1] = 0.147, lam_tx[2] = 0.45; //k_1, k_2
	lam_tx[3] = 2.465, lam_tx[4] = 1.806, lam_tx[5] = 0.906;  //a_1, a_2, a_3 
	lam_tx[6] = 0.431, lam_tx[7] = 0.39, lam_tx[8] = 0.0626; //z0_1, z0_2, z_03
	lam_tx[9] = 0.0913, lam_tx[10] = 0.636, lam_tx[11] = 0.123; //k_m1, k_m2, k_m3
	
	//Set histo to 0...................................................................
	for(i=0; i<bins_m; i++) for (j=0; j<bins_z; j++) for (k=0; k<bins_t; k++) histo[i][j][k] = 0;

	//Buscar extrems de la mostra......................................................
	mock = file_init(cat_file);

	N = 0;
	fscanf(mock, "%f %f %f\n", &m, &z, &t);

	m_max = m;
	zmax = z;

	m_min = m;
	zmin = z;

	N++;
	while(fscanf(mock, "%f %f %f\n", &m, &z, &t) != -1){
		if(m > m_max) m_max=m;
		if(z > zmax) zmax=z;
		
		if(m < m_min) m_min=m;
		if(z < zmin) zmin=z;
		
		N++; }
		
	printf("\n---PLOTTING PRIOR MODULE------------------------------------\n");
	printf("\tTotal num. gal. in cat.:\t%d\n", N);
	printf("\tTotal Magnitude range:\t\t[%2.2f,%2.2f]\n", m_min, m_max);
	printf("\tTotal Redshift range:\t\t[%2.2f,%2.2f]\n", zmin, zmax);
	printf("\tMagnitude bins:\t\t\t%d\n", bins_m);
	printf("\tRedshift bins:\t\t\t%d\n", bins_z);
	printf("\tType bins:\t\t\t%d\n", bins_t);

	fclose(mock);

	//Construir histograma................................................................
	mock=file_init(cat_file);
	l = 0, n_fora=0;

	delta_m=m_max-m_min;
	delta_z=zmax-zmin;

	while(fscanf(mock, "%f %f %f\n", &m, &z, &t) != -1){
		l++;
			
		if(t <= t1) t = 1;
		else{
			if(t <= t2) t = 2;
			else t = 3;}
		
		i=(int)(bins_m*((m-m_min)/delta_m));
		j=(int)(bins_z*((z-zmin)/delta_z));
		k=t-1;
		
		if(i < bins_m && i >= 0 && j < bins_z && j >= 0 && k < bins_t && k>=0 ) histo[i][j][k]++;
		else n_fora++;}
		
		//printf("i=%d j=%d k=%d\n", i, j, k);
		//printf("l=%d | m=%f z=%f t=%f |",l ,m ,z ,t);
		//printf(" histo[%d][%d][%d]=%f\n", i, j, k, histo[i][j][k]);


	//Check de N...........................................................................
	l = 0;
	for(i=0; i<bins_m; i++) for(j=0; j<bins_z; j++) for(k=0; k<bins_t; k++) l+=(int)histo[i][j][k];	
	if(l+n_fora != N){ 
		printf("\n\tError: Total num. of gal. don't match with num. in histo (%d+%d!=%d)\n",l,n_fora,N);
		exit(0);}

	//Normalització.......................................................................
	for(i=0; i<bins_m; i++) for (j=0; j<bins_z; j++) for (k=0; k<bins_t; k++) histo[i][j][k]=(float)(histo[i][j][k]/(l*(delta_m/bins_m)*(delta_z/bins_z)));

	//Check de Normalització..............................................................
	norm = 0;
	for(i = 0; i<bins_m; i++) for(j=0; j<bins_z; j++)  for(k=0; k<bins_t; k++) norm+=histo[i][j][k]*(delta_m/bins_m)*(delta_z/bins_z);
	if(fabs(norm-1)>0.0001){ 
		printf("\n\tError: Normalization hasn't been computed correctly (%f!=1.0)\n",norm);
		exit(0);}
	
	printf("\n\tCreating plots and converting to pdf...\n");
	for(m_prior = m_0_tx; m_prior <= m_cut; m_prior += 0.5){ //Loop over magnitude
		dm_tx = m_prior - m_0_tx, dm = m_prior - m_0;
	
		//Fer el plot................................................................
		plot = fopen("plot.txt", "w");
		fprintf(plot, "set terminal postscript enhanced color\nset output \"%sp_z_t_%1.1f.ps\"\nset encoding iso_8859_1\nset multiplot\nset size 0.5,0.5\n", plot_file, m_prior);

		//Compute p(t|m)..............................................................
		sprintf(nom_fich_tm, "%sp_t_%1.1f.txt", plot_txt_file, m_prior);
		sprintf(nom_fich_tm_func, "%sp_t_%1.1f_func.txt", plot_txt_file, m_prior);
		sprintf(nom_fich_tm_func_tx, "%sp_t_%1.1f_func_tx.txt", plot_txt_file, m_prior);

		p_t_m_file=fopen(nom_fich_tm, "w");
		p_t_m_func_file=fopen(nom_fich_tm_func, "w");
		p_t_m_func_tx_file=fopen(nom_fich_tm_func_tx, "w");

		i=(int)(bins_m*((m_prior-m_min)/delta_m));
		
		//Normalization.....................................................
		p_tm[0]=0;
		for(k=0; k<bins_t; k++) for(j=0; j<bins_z; j++) p_tm[0]+=histo[i][j][k]*(delta_z/bins_z);
		
		//Writing files with curves..........................................
		for (k=0; k<bins_t; k++){
			p_tm[1]=0;

			for (j=0; j<bins_z; j++) p_tm[1]+=histo[i][j][k]*(delta_z/bins_z);
			
			fprintf(p_t_m_file,"%d %f\n", k+1 , p_tm[1]/p_tm[0]);
			if(k+1 < 3){
				fprintf(p_t_m_func_file,"%d %f\n", k+1 , f[k+1]*exp(-fabs(lam[k+1])*dm));
				fprintf(p_t_m_func_tx_file,"%d %f\n", k+1 , f_tx[k+1]*exp(-fabs(lam_tx[k+1])*dm_tx));}
			else{
				fprintf(p_t_m_func_file,"%d %f\n", k+1 , 1-f[1]*exp(-fabs(lam[1])*dm)-f[2]*exp(-fabs(lam[2])*dm));
				fprintf(p_t_m_func_tx_file,"%d %f\n", k+1 , 1-f_tx[1]*exp(-fabs(lam_tx[1])*dm_tx)-f_tx[2]*exp(-fabs(lam_tx[2])*dm_tx));
			}
		}
			
		fclose(p_t_m_file), fclose(p_t_m_func_file), fclose(p_t_m_func_tx_file);

		fprintf(plot, "set xlabel \"{t}\"\nset ylabel \"{p(t|%1.1f)}\"\nset origin 0.0, 0.5\nplot [0.5:3.5] [0:1] \"%s\" title \"{/Helvetica=10 Real values}\" with histeps linewidth 3, \"%s\" title \"{/Helvetica=10 Fit values}\" with linespoints lw 3, \"%s\" title \"{/Helvetica=10 Txitxo values}\" with linespoints lw 3\n", m_prior,nom_fich_tm, nom_fich_tm_func, nom_fich_tm_func_tx);
	
		//Calcul p(z|t,m)........................................................
		origin[0][0] = 0.5, origin[0][1] = 0.5;
		origin[1][0] = 0.0, origin[1][1] = 0.0;
		origin[2][0] = 0.5, origin[2][1] = 0.0;

		for(t_prior = 1; t_prior <= 3; t_prior++){ //Loop over type
			sprintf(nom_fich_ztm, "%sp_z_%d_%1.1f.txt", plot_txt_file, t_prior, m_prior);
			sprintf(nom_fich_ztm_func, "%sp_z_%d_%1.1f_func.txt", plot_txt_file, t_prior, m_prior);
			sprintf(nom_fich_ztm_func_tx, "%sp_z_%d_%1.1f_func_tx.txt", plot_txt_file, t_prior, m_prior);
			
			p_z_t_m_file = fopen(nom_fich_ztm, "w");
			p_z_t_m_func_file = fopen(nom_fich_ztm_func, "w");
			p_z_t_m_func_tx_file = fopen(nom_fich_ztm_func_tx, "w");

			//Computing p(z|t,m)................................................
			zmt = fabs(lam[t_prior+5]) + (fabs(lam[t_prior+8])*dm);
			zmt_tx = fabs(lam_tx[t_prior+5]) + (fabs(lam_tx[t_prior+8])*dm_tx);

			//Normalization.....................................................
			zi = 0, sum = 0, dz = 0.01; //Fitted
			while(zi<=20){
				sum += (dz*pow(zi,fabs(lam[t_prior+2]))*exp(-pow(zi/zmt,fabs(lam[t_prior+2]))));
				zi += dz;}

			zi = 0, sum_tx = 0; //Txitxo
			while(zi <= 20){
				sum_tx += (dz*pow(zi,fabs(lam_tx[t_prior+2]))*exp(-pow(zi/zmt_tx,fabs(lam_tx[t_prior+2]))));
				zi += dz;}
			
			p_ztm = 0, i = (int)(bins_m*((m_prior-m_min)/delta_m)), k = t_prior-1; //Current
			for (j=0; j<bins_z; j++) p_ztm+=histo[i][j][k]*(delta_z/bins_z);
			
			//Writing files with curves..........................................
			for (j=0; j<bins_z; j++){
				z = (float)( (((j+0.5)*delta_z)/bins_z)+zmin );	
				fprintf(p_t_m_file,"%f %f\n", z, histo[i][j][k]/p_ztm);
				fprintf(p_t_m_func_file,"%f %f\n", z , (1/sum)*(pow(z,fabs(lam[t_prior+2]))*exp(-pow(z/zmt,fabs(lam[t_prior+2])))) );
				fprintf(p_t_m_func_tx_file,"%f %f\n", z , (1/sum_tx)*(pow(z,fabs(lam_tx[t_prior+2]))*exp(-pow(z/zmt_tx,fabs(lam_tx[t_prior+2])))) );}
			
			fclose(p_z_t_m_file), fclose(p_z_t_m_func_file), fclose(p_z_t_m_func_tx_file);
			
			fprintf(plot, "set xlabel \"{z}\"\nset xrange[0:4]\nset ylabel \"{p(z|%d, %1.1f)}\"\nset origin %f, %f\nplot \"%s\" title \"{/Helvetica=10 Real values}\" with histeps linewidth 3, \"%s\" title \"{/Helvetica=10 Fit values}\" with linespoints lw 3, \"%s\" title \"{/Helvetica=10 Txitxo values}\" with linespoints lw 3\n", t_prior, m_prior, origin[t_prior-1][0], origin[t_prior-1][1], nom_fich_ztm, nom_fich_ztm_func, nom_fich_ztm_func_tx);
			//fprintf(plot, "!pstopdf %sp_z_t_%1.1f.ps %sp_z_t_%1.1f.pdf", plot_file, m_prior, plot_file, m_prior);
		}
		fclose(plot);

		system("gnuplot \"plot.txt\" ");
		system("rm \"plot.txt\" ");
	}
}

void powell(float p[], float **xi, int n, float ftol, int *iter, float *fret, float (*func)(float [])) //Minimization of a function func of n variables. Input consists of an initial starting point p[1..n]; an initial matrix xi[1..n][1..n], whose columns contain the initial set of directions (usually the n unit vectors); and ftol, the fractional tolerance in the function value such that failure to decrease by more than this amount on one iteration signals doneness. On output, p is set to the best point found, xi is the then-current direction set, fret is the returned function value at p, and iter is the number of iterations taken. The routine linmin is used.
{
	void linmin(float p[], float xi[], int n, float *fret, float (*func)(float []));
	int i,ibig,j;
	float del,fp,fptt,t,*pt,*ptt,*xit;
	pt=vector(1,n);
	ptt=vector(1,n);
	xit=vector(1,n);
	*fret=(*func)(p);
	for (j=1;j<=n;j++) pt[j]=p[j]; //Save the initial point
	for (*iter=1;;++(*iter)) 
	{
		fp=(*fret);
		ibig=0;
		del=0.0; //Will be the biggest function decrease
		for (i=1;i<=n;i++) //In each iteration, loop over all directions in the set
		{ 
			for (j=1;j<=n;j++) 
			{
				xit[j]=xi[j][i]; //Copy the direction
			}
			
			fptt=(*fret);
			linmin(p,xit,n,fret,func); //minimize along it
			if (fabs(fptt-(*fret)) > del) //and record it if it is the largest decrease so far.
			{ 	
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) 
		{
			free_vector(xit,1,n); //Termination criterion.
			free_vector(ptt,1,n);
			free_vector(pt,1,n);
			return;
		}
		if (*iter == ITMAX) nrerror("powell exceeding maximum iterations.");
		for (j=1;j<=n;j++) //Construct the extrapolated point and the average direction moved. Save the old starting point.
		{ 
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt); //Function value at extrapolated point.
		if (fptt < fp) 
		{
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0) 
			{
				linmin(p,xit,n,fret,func); //Move to the minimum of the new direction, and save the new direction.
				for (j=1;j<=n;j++) 
				{
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	} //Back for another iteration.
}

int ncom; //Global variables communicate with f1dim.
float *pcom,*xicom,(*nrfunc)(float []);

void linmin(float p[], float xi[], int n, float *fret, float (*func)(float [])) //Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and resets p to where the function func(p) takes on a minimum along the direction xi from p, and replaces xi by the actual vector displacement that p was moved. Also returns as fret the value of func at the returned location p. This is actually all accomplished by calling the routines mnbrak and brent.
{
	float brent(float ax, float bx, float cx, float (*f)(float), float tol, float *xmin);
	float f1dim(float x);
	void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc, float (*func)(float));
	int j;
	float xx,xmin,fx,fb,fa,bx,ax;
	ncom=n; //Dene the global variables.
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) 
	{
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0; //Initial guess for brackets.
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=1;j<=n;j++) //Construct the vector results to return. 
	{ 
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	
	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
}

extern int ncom; //Dened in linmin.
extern float *pcom,*xicom,(*nrfunc)(float []);

float f1dim(float x) //Must accompany linmin.
{
	int j;
	float f,*xt;
	xt=vector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_vector(xt,1,ncom);
	return f;
}


#define NRANSI
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc, float (*func)(float))
{
	float ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
	  r=(*bx-*ax)*(*fb-*fc);
	  q=(*bx-*cx)*(*fb-*fa);
	  u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
	    (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
	  ulim=(*bx)+GLIMIT*(*cx-*bx);
	  if ((*bx-u)*(u-*cx) > 0.0) {
	    fu=(*func)(u);
	    if (fu < *fc) {
	      *ax=(*bx);
	      *bx=u;
	      *fa=(*fb);
	      *fb=fu;
	      return;
	    } else if (fu > *fb) {
	      *cx=u;
	      *fc=fu;
	      return;
	    }
	    u=(*cx)+GOLD*(*cx-*bx);
	    fu=(*func)(u);
	  } else if ((*cx-u)*(u-ulim) > 0.0) {
	    fu=(*func)(u);
	    if (fu < *fc) {
	      SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
		SHFT(*fb,*fc,fu,(*func)(u))
		}
	  } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
	    u=ulim;
	    fu=(*func)(u);
	  } else {
	    u=(*cx)+GOLD*(*cx-*bx);
	    fu=(*func)(u);
	  }
	  SHFT(*ax,*bx,*cx,u)
	  SHFT(*fa,*fb,*fc,fu)
       }
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef NRANSI

//#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
//Here ITMAX is the maximum allowed number of iterations; CGOLD is the golden ratio; ZEPS is a small number that protects against trying to achieve fractional accuracy for a minimum that happens to be exactly zero.
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

float brent(float ax, float bx, float cx, float (*f)(float), float tol, float *xmin) //Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates the minimum to a fractional precision of about tol using Brent's method. The abscissa of the minimum is returned as xmin, and the minimum function value is returned as brent, the returned function value.
{
	int iter;
	float a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	float e=0.0; //This will be the distance moved on the step before last.
	a=(ax < cx ? ax : cx); //a and b must be in ascending order,
	b=(ax > cx ? ax : cx); //but input abscissas need not be.
	x=w=v=bx; //Initializations...
	fw=fv=fx=(*f)(x);
	for (iter=1;iter<=ITMAX;iter++) //Main program loop. 
	{ 
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) //Test for done here.
		{			
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) //Construct a trial parabolic t. 
		{
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) 
			{
				d=CGOLD*(e=(x >= xm ? a-x : b-x)); //The above conditions determine the acceptability of the parabolic fit. Here we take the golden section step into the larger of the two segments.
			}
			else 
			{
				d=p/q; //Take the parabolic step.
				u=x+d;
				if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
			}
		} 
		else 
		{
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		//This is the one function evaluation per iteration.
		if (fu <= fx) //Now decide what to do with our function evaluation. 
		{ 
			if (u >= x) a=x; else b=x; 
			SHFT(v,w,x,u) //Housekeeping follows:
			SHFT(fv,fw,fx,fu)
		} 
		else 
		{
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) 
			{
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} 
			else if (fu <= fv || v == x || v == w) 
			{
				v=u;
				fv=fu;
			}
		} //Done with housekeeping. Back for another iteration.
	}	
	nrerror("Too many iterations in brent");
}



