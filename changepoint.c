#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "pcg_variants.h"
#include "gauss.h"

#include <time.h>
#include <stdlib.h>


#define RAND_MAX_U32 4294967295
#define MOVE   (uint32_t) (0.20 * (RAND_MAX_U32+1))
#define UPDATE (uint32_t) (0.15 * (RAND_MAX_U32+1))
#define SIGMA  (uint32_t) (0.15 * (RAND_MAX_U32+1))
#define BIRTH  (uint32_t) (0.25 * (RAND_MAX_U32+1))
#define DEATH  (uint32_t) (0.25 * (RAND_MAX_U32+1))



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
	const uint32_t K=rows-1; // Number of possible changepoint locations
	const uint32_t npmin=0; // Minimum number of changepoints
	const uint32_t npmax=K; // Maximum number of changepoints
	const uint32_t dmin=10; // Minimum number of points needed between any two changepoints (in any partition)

	const uint32_t nsims = (uint32_t)1E9; // Number of simulations to run
	printf("Nsims: %lu\n",nsims);
	
	// De-mean, normalize variance, and define initial model
	for(j=0;j<columns;j++){
		standardize(&d[j*rows], rows);
		for(i=0; i<rows; i++){
			s[j*rows+i]=1.0;
			m[j*rows+i]=0.0;
		}
	}

	// Initialize random number generator
	pcg32_random_t rng;
	pcg32_srandom_r(&rng,time(NULL), clock());
	uint32_t r;
	double a;

	clock_t t0, t1, t2, t3, t4;


	// Make arrays to hold number of changepoints, r-squared, and log likelihood
	uint32_t * const restrict np=malloc(sizeof(int)*nsims);
	double * const restrict r2C=malloc(sizeof(double)*nsims);
	double * const restrict llC=malloc(sizeof(double)*nsims);
	np[0]=0;

	// Make array to hold boundary points between partitions
	uint32_t * const restrict bndPnts=malloc(sizeof(uint32_t)*(K+2));
	bndPnts[0]=0;
	bndPnts[np[0]+1]=rows;
	for (i=0; i<np[0]; i++){
		bndPnts[i+1]=pcg32_random_r(&rng)/(RAND_MAX_U32/K);
	}
	np[0]=unique_uints(bndPnts,np[0]+2)-2;


	// The actual loop
	for (i=0;i<nsims;i++){

		// Randomly choose a a modification to the model 
		r = pcg32_random_r(&rng);
		a = pcg32_random_r(&rng)/4294967295.0;
//		rn = pcg_gaussian_ziggurat(&rng, 1.0);

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

	return 0;
}

double likelihood(double* const bndPnts){
	return 0;
}
