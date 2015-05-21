#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "arrays.h"
#include "pcg_variants.h"
#include "gauss.h"

#define RAND_MAX_U32 4294967295
#define MOVE   (uint32_t) (0.30 * (RAND_MAX_U32+1))
#define UPDATE (uint32_t) (0.00 * (RAND_MAX_U32+1))
#define SIGMA  (uint32_t) (0.20 * (RAND_MAX_U32+1))
#define BIRTH  (uint32_t) (0.25 * (RAND_MAX_U32+1))
#define DEATH  (uint32_t) (0.25 * (RAND_MAX_U32+1))

double log_likelihood(const double* const d, const double* const m, const double* const s, const uint32_t rows, const uint32_t columns){
	double ll=0;
	for (uint32_t i=0; i<columns; i++){
		for (uint32_t j=0; j<rows; j++){
			ll -= pow((d[i*rows+j]-m[i*rows+j])/s[i],2);
		}
	}
	return ll/2;
}

void update_model(const uint32_t* const bndPntsC, const uint32_t np, const uint32_t rows, const uint32_t columns, pcg32_random_t* rng, const double* const d, const double* const s, double* const restrict m){
	uint32_t i, j, k;
	double mu, sigma, A;
	for (i=0; i<columns; i++){
//		min=minArray(&d[i*rows],rows);
//		max=maxArray(&d[i*rows],rows);
		for (j=0; j<(np+1); j++){
//			min=minArray(&d[i*rows + bndPntsC[j]],bndPntsC[j+1]-bndPntsC[j]);
//			max=maxArray(&d[i*rows + bndPntsC[j]],bndPntsC[j+1]-bndPntsC[j]);
			Offset_nanstd(&d[i*rows + bndPntsC[j]], bndPntsC[j+1]-bndPntsC[j], &mu, &sigma);
//			A = (double) pcg32_random_r(rng) / 42949672959.0 * (max-min) + min;
			A = pcg_gaussian_ziggurat(rng, sigma) + mu;
			for (k=bndPntsC[j]; k<bndPntsC[j+1]; k++){
				m[i*rows+k]=A;
			}
		}
	}
}

void update_sigma(const uint32_t* const bndPntsC, const uint32_t np, const uint32_t rows, const uint32_t columns, const double* const d, double* const restrict m){

}


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
	double* const restrict s = malloc(sizeof(double)*columns);
	double* const restrict m = malloc(sizeof(double)*rows*columns);

	double mu[columns], sigma[columns];

	// Markov chain version using data uncertainties
	const uint32_t K=rows-1; // Number of possible changepoint locations
	const uint32_t npmin=0; // Minimum number of changepoints
	const uint32_t npmax=K; // Maximum number of changepoints
	const uint32_t dmin=10; // Minimum number of points needed between any two changepoints (in any partition)

	const uint32_t nsims = 1E5; // Number of simulations to run
	printf("Nsims: %u\n",nsims);
	
	// De-mean, normalize variance, and define initial model
	for(j=0;j<columns;j++){
		standardize(&d[j*rows], rows);
		s[j]=1;
		for(i=0; i<rows; i++){
			m[j*rows+i]=0.0;
		}
	}

	// Initialize random number generator
	pcg32_random_t rng;
	pcg32_srandom_r(&rng,time(NULL), clock());
	uint32_t r;
	double u;

	// Make arrays to hold number of changepoints, r-squared, and log likelihood
//	uint32_t * const restrict np=malloc(sizeof(int)*nsims);
	uint32_t np, npP;
//	double * const restrict r2C=malloc(sizeof(double)*nsims);
//	double * const restrict llC=malloc(sizeof(double)*nsims);
	np=2;

	// Make array to hold boundary points between partitions
	uint32_t * const restrict bndPntsC=malloc(sizeof(uint32_t)*(K+2));
	uint32_t * const restrict bndPntsP=malloc(sizeof(uint32_t)*(K+2));
	bndPntsC[0]=0;
	bndPntsC[np+1]=rows;
	for (i=0; i<np; i++){
		bndPntsC[i+1]=pcg32_random_r(&rng)/(RAND_MAX_U32/K);
	}
	np=unique_uints(bndPntsC,np+2)-2;
	copyArrayUint(bndPntsC,K+2,bndPntsP);


	double llP, ll=log_likelihood(d, m, s, rows, columns);


	// The actual loop
	for (i=0;i<nsims;i++){

		// Randomly choose a a modification to the model 
		r = pcg32_random_r(&rng);
		u = log(pcg32_random_r(&rng)/4294967295.0);
		uint32_t pick;	

		if (r < MOVE){
			if (np>0){
				copyArrayUint(bndPntsC,K+2,bndPntsP);
				// Pick which changepoint to move
				pick=pcg32_random_r(&rng)/(RAND_MAX_U32/np)+1;
				// Move the changepoint between its boundaries
				bndPntsP[pick] = pcg32_random_r(&rng) / ( RAND_MAX_U32 / ((bndPntsP[pick+1]-1)-(bndPntsP[pick-1]+1)+1) ) + (bndPntsP[pick-1]+1);

				// Update the model
				update_model(bndPntsP, np, rows, columns, &rng, d, s, m);

				// Calculate log likelihood for proposal
				llP=log_likelihood(d, m, s, rows, columns);
				if (u<llP-ll){
					ll=llP;
					copyArrayUint(bndPntsP,K+2,bndPntsC);
				}
			}
		} else if (r < MOVE+UPDATE){
			// Update the model
			update_model(bndPntsP, np, rows, columns, &rng, d, s, m);

			// Calculate log likelihood for proposal
			llP=log_likelihood(d, m, s, rows, columns);
			if (u<llP-ll){
				ll=llP;
				copyArrayUint(bndPntsP,K+2,bndPntsC);
			}

		} else if (r < MOVE+UPDATE+SIGMA){
			// Choose a new sigma
			for(j=0; j<columns; j++){
				s[j]=log(pcg32_random_r(&rng)/4294967295.0);
			}
			
			// Calculate log likelihood for proposal
			llP=log_likelihood(d, m, s, rows, columns);
			if (u<llP-ll){ // If accepted
				ll=llP;
				copyArrayUint(bndPntsP,K+2,bndPntsC);
			}

		} else if (r < MOVE+UPDATE+SIGMA+BIRTH){
			if (np < npmax) {
				copyArrayUint(bndPntsC,K+2,bndPntsP);

				// Pick which changepoint to add right of
				pick=pcg32_random_r(&rng)/(RAND_MAX_U32/(np+1));
				bndPntsP[np+2] = pcg32_random_r(&rng) / ( RAND_MAX_U32 / ((bndPntsP[pick+1]-1)-(bndPntsP[pick]+1)+1) ) + (bndPntsP[pick]+1);
				npP=unique_uints(bndPntsP,np+3)-2;

				// Update the model
				update_model(bndPntsP, npP, rows, columns, &rng, d, s, m);
				// Calculate log likelihood for proposal
				llP=log_likelihood(d, m, s, rows, columns);
				if (u<llP-ll){
					ll=llP;
					np=npP;
					copyArrayUint(bndPntsP,K+2,bndPntsC);
				}
			}
		} else if (r < MOVE+UPDATE+SIGMA+BIRTH+DEATH){
			if (np > npmin) {
				copyArrayUint(bndPntsC,K+2,bndPntsP);


				// Pick which changepoint to delete
				pick=pcg32_random_r(&rng)/(RAND_MAX_U32/np)+1;
				bndPntsP[pick]=bndPntsP[pick+1];
				npP=unique_uints(bndPntsP,np+2)-2;

				// Update the model
				update_model(bndPntsP, np, rows, columns, &rng, d, s, m);
				// Calculate log likelihood for proposal
				llP=log_likelihood(d, m, s, rows, columns);
				if (u<llP-ll){
					ll=llP;
					np=npP;
					copyArrayUint(bndPntsP,K+2,bndPntsC);
				}
			}
		} else {
			printf("Uh-oh!\n");
		}
//		printf("%u\n",i);
	}

	return 0;
}


