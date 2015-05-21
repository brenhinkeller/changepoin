#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "mt.h"

#include <time.h>
#include <stdlib.h>


#define MT_RAND_MAX 4294967295
#define MOVE   (uint32_t) (0.2 * MT_RAND_MAX)
#define UPDATE (uint32_t) (0.15 * MT_RAND_MAX)
#define SIGMA  (uint32_t) (0.15 * MT_RAND_MAX)
#define BIRTH  (uint32_t) (0.25 * MT_RAND_MAX)
#define DEATH  (uint32_t) (0.25 * MT_RAND_MAX)



int main(int argc, char **argv){
	// Check input arguments
	if (argc != 2) {
		printf("USAGE: %s <input_filename>\n", argv[0]);
		exit(1);
	}

	uint32_t rows, columns;
	uint32_t i, j;
	// import 2-d data array as a flat double array
	double* const restrict d = csvparseflat(argv[1],',', &rows, &columns);		
	double* const restrict s = malloc(sizeof(double)*rows*columns);
	double* const restrict m = malloc(sizeof(double)*rows*columns);

	double mu[columns], sigma[columns];

	// Markov chain version using data uncertainties
	const uint32_t K=100;//=rows-1; // Number of possible changepoint locations
	const uint32_t npmin=0; // Minimum number of changepoints
	const uint32_t npmax=K; // Maximum number of changepoints
	const uint32_t dmin=10; // Minimum number of points needed between any two changepoints (in any partition)

	const uint32_t nsims = 10000; // Number of simulations to run


	// De-mean, normalize variance, and define initial model
	for(j=0;j<columns;j++){
		standardize(&d[j*rows], rows);
		for(i=0; i<rows; i++){
			s[j*rows+i]=1.0;
			m[j*rows+i]=0.0;
		}
	}

	// Initialize random number generator
	uint32_t mt_buffer[MT_LEN], mt_index;
	double r;
//	int mt_index;
	srandom(time(NULL));
	mt_init(mt_buffer,&mt_index);


	clock_t t0, t1, t2, t3, t4;
	

	// Make arrays to hold number of changepoints, r-squared, and log likelihood
	uint32_t * const restrict np=malloc(sizeof(int)*nsims);
	np[0]=100;
	double * const restrict r2C=malloc(sizeof(double)*nsims);

	double * const restrict llC=malloc(sizeof(double)*nsims);


	// Make array to hold boundary points between partitions
	uint32_t * const restrict bndPnts=malloc(sizeof(uint32_t)*(K+2+1000));
	bndPnts[0]=0;
	for (i=0; i<np[0]; i++){
		bndPnts[i+1]=mt_rand_r(mt_buffer,&mt_index)/(MT_RAND_MAX/100);
		printf("%g,",mt_rand_r(mt_buffer,&mt_index)/(double) MT_RAND_MAX);
	}
	printf("\n");
//	unique_uints(&bndPnts[1],&np[0]);
	bndPnts[np[0]+2]=rows;

	printf("Number of changepoints: %u\n", np[0]);
	for (i=0; i<np[0]+1; i++){
		printf("%i,",bndPnts[i+1]);
	}
	printf("\n");

	// The actual loop
	for (i=0;i<nsims;i++){

		// Randomly choose a a modification to the 
		r = mt_rand_r(mt_buffer, &mt_index);

		if (r < MOVE){


		} else if (r < MOVE+UPDATE){


		} else if (r < MOVE+UPDATE+SIGMA){


		} else if (r < MOVE+UPDATE+SIGMA+BIRTH){

			if (np[i] < npmax) {
			}
		} else if (r < MOVE+UPDATE+SIGMA+BIRTH+DEATH){

			if (np[i] > npmin) {
			}
		} else {
			printf("Uh-oh!\n");
		}

	}


	//Calculate mean and variance
	double mean[columns], var[columns];
	t1=clock();
	for(j=0;j<columns;j++){
		Knuth_nanstd(&d[j*rows], rows, &mean[j], &var[j]);
	}
	t2=clock();

	//Calculate mean and variance
	for(j=0;j<columns;j++){
		var[j]=nanstd(&d[j*rows], rows);
	}	
	t3=clock();

	//Calculate mean and variance
	for(j=0;j<columns;j++){
		mean[j]=nanmean(&d[j*rows], rows);
	}	
	t4=clock();


	printf("Elapsed time (Knuth_nanstd): %lu\n",t2-t1);
	printf("Elapsed time (nanstd): %lu\n",t3-t2);
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

double likelihood(double* const bndPnts){
	return 0;
}
