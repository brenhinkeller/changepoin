#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "mt.h"

#include <time.h>
#include <stdlib.h>


int main(int argc, char **argv){
	// Check input arguments
	if (argc != 2) {
		printf("USAGE: %s <input_filename>\n", argv[0]);
		exit(1);
	}

	int rows, columns;
	int i, j;
	// import 2-d data array as a flat double array
	double* const restrict data=csvparseflat(argv[1],',', &rows, &columns);

		
	double* d = malloc(sizeof(double)*rows*columns);
	double* s = malloc(sizeof(double)*rows*columns);
	double* m = malloc(sizeof(double)*rows*columns);

	double mu[columns], sigma[columns];

	// Markov chain version using data uncertainties
	//
	const int possiblepoints=rows-1; // Number of possible changepoint locations
	const int npmin=0; // Minimum number of changepoints
	const int npmax=possiblepoints; // Maximum number of changepoints
	const int dmin=10; // Minimum number of points needed between any two changepoints (in any partition)

	const int nsims = 10000; // Number of simulations to run



	//Calculate mean and variance
	for(j=0;j<columns;j++){
		nanstd(&data[j*rows], rows, &mu[j], &sigma[j]);
	}

	// De-mean, normalize variance, and define initial model
	for(j=0;j<columns;j++){
		for(i=0; i<rows; i++){
			d[j*rows+i]=(data[j*rows+i]-mu[j])/sigma[j];
			s[j*rows+i]=1.0;
			m[j*rows+i]=0.0;
		}
	}


	clock_t t0, t1, t2, t3, t4;

	printf("%g\n",0.0/0);
	

	// Propose a random number of changepoints
	// np=round(rand(1,1)*(npmax-npmin)+npmin);
	int *np=malloc(sizeof(int)*nsims);
	np[0]=0;
	double *r2C=malloc(sizeof(double)*nsims);
//	r2C(1)=NaN; % Current sum squared residual
	double *llC=malloc(sizeof(double)*nsims);
//	llC[0]=

	int *bndPnts=malloc(sizeof(int)*(possiblepoints+2));

	sort_ints(bndPnts,possiblepoints+2);

//	% Propose an initial random set of boundary points
//	bndPntsC=[0; sort(ceil(rand(np(1),1)*possiblepoints)); size(data,1)];
//
//	% means=NaN(npmax+2,1);
//	meansP=NaN(length(bndPntsC)-1,1);
//

	printf("%i\n", nsims);



//	t1=clock();
//		unsigned long mt_buffer[MT_LEN];
//		int mt_index;
//		mt_init(mt_buffer,&mt_index);	
//		double* r1=malloc(1000000000*sizeof(double));
//		for (i=0; i<1000000000; i++){r1[i]=mt_rand_r(mt_buffer, &mt_index)/4294967295.0;}
//	t2=clock();
//
//	printf("Elapsed time (rand_r): %lu\n",t2-t1);


//	// Print the array
//	for (i=0;i<rows;i++){
//		for (j=0;j<columns;j++){
//			printf("%g, ",data[j*rows + i]);
//		}
//		printf("\n");
//	}
	


	//Calculate mean and variance
	double mean[columns], var[columns];
	t1=clock();
	for(j=0;j<columns;j++){
		Knuth_var(&data[j*rows], rows, &mean[j], &var[j]);
	}
	t2=clock();

	//Calculate mean and variance
	for(j=0;j<columns;j++){
		nanvar(&data[j*rows], rows, &mean[j], &var[j]);
	}	
	t3=clock();

	//Calculate mean and variance
	for(j=0;j<columns;j++){
		mean[j]=nanmean(&data[j*rows], rows);
	}	
	t4=clock();


	printf("Elapsed time (Knuth var): %lu\n",t2-t1);
	printf("Elapsed time (nanvar): %lu\n",t3-t2);
	printf("Elapsed time (nanmean): %lu\n",t4-t3);

	//Print mean and variance
	for(j=0;j<columns;j++){
		printf("%g,",mean[j]);
	}
	printf("\n");
	for(j=0;j<columns;j++){
		printf("%g,",var[j]);
	}
	printf("\n");


	normalize(mean,columns);
	for(j=0;j<columns;j++){
		printf("%g,",mean[j]);
	}
	printf("\n");


	return 0;
}


