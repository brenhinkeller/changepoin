#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

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


double lqxz(const double * const d, const double * const m, const uint32_t * const partition, const uint32_t rows, const uint32_t columns){
	double mu, sigma;
	double lqz=0, lq=0, lz=0;
	bool exists=0;
	for (uint32_t i=0; i<columns; i++){
		Offset_nanstd(&d[i*rows + partition[0]], partition[1]-partition[0], &mu, &sigma);
		sigma+=1E-15;
		lq = -log(sigma*sqrt(2*3.1415926535)) - pow( (d[i*rows + partition[0]] - mu)/(sigma) , 2) / 2;
		lz = log(maxArray(&d[i*rows],rows) - minArray(&d[i*rows],rows) + 1E-15); 
		if (!isnan(lq) && !isnan(lz)){
			lqz += lq + lz;
			exists=1;
		}
//		printf("mu: %g, sigma: %g\n", mu, sigma);
//		printf("lq: %g, lz: %g\n", lq, lz);
	}
//	printf("partition[0]: %u, partition[1]: %u\n", partition[0], partition[1]);
	if (exists){
		return lqz;
	} else {
		return NAN;
	}
}




// Currently not used
//void update_sigma(const uint32_t* const bndPntsC, const uint32_t np, const uint32_t rows, const uint32_t columns, const double* const d, double* const restrict m){
//}


int main(int argc, char **argv){
	// Check input arguments
	if (argc != 3) {
		fprintf(stderr,"USAGE: %s <number of simulations> <input_filename>\n", argv[0]);
		exit(1);
	}

	uint32_t rows, columns;
	uint32_t i, j;
	// import 2-d data array as a flat double array
	double* const restrict d = csvparseflat(argv[2],',', &rows, &columns);		
	double* const restrict s = malloc(sizeof(double)*columns);
	double* const restrict sP = malloc(sizeof(double)*columns);
	double* const restrict m = malloc(sizeof(double)*rows*columns);

	double mu[columns], sigma[columns];

	// Markov chain version using data uncertainties
	const uint32_t K=rows-1; // Number of possible changepoint locations
	const uint32_t npmin=0; // Minimum number of changepoints
	const uint32_t npmax=K; // Maximum number of changepoints
	const uint32_t dmin=10; // Minimum number of points needed between any two changepoints (in any partition)

	const uint32_t nsims = atoi(argv[1]); // Number of simulations to run
	fprintf(stderr, "Nsims: %u\n",nsims);
	
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
//	double * const restrict r2C=malloc(sizeof(double)*nsims);
//	double * const restrict llC=malloc(sizeof(double)*nsims);
	uint32_t np=2, npP;

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


	uint32_t pick;	
	double lqz, sll, llP, ll=log_likelihood(d, m, s, rows, columns);

	// The actual loop
	for (i=0;i<nsims;i++){


		// Randomly choose a type of modification to the model 
		r = pcg32_random_r(&rng);
		u = log(pcg32_random_r(&rng)/4294967295.0);


		// Move a changepoint
		if (r < MOVE && np>0){ 
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
//				printf("Move: llP-ll = %g - %g\n",llP,ll);
//				printf("Accepted!\n");
				ll=llP;
				copyArrayUint(bndPntsP,K+2,bndPntsC);
//				for (int n=0; n<np+2; n++){
//					printf("%u,",bndPntsP[n]);
//				}
//				printf("\n%u\n",np);
			}
		} 


		// Adjust SIGMA
		else if (r < MOVE+SIGMA){
			copyArray(s,columns,sP);

			// Choose a new sigma
			j=pcg32_random_r(&rng)/(RAND_MAX_U32/columns);
			sP[j]=pcg32_random_r(&rng)/4294967295.0;
			sll=log(pow((sP[j]/s[j]),rows));

			// Calculate log likelihood for proposal
			llP=log_likelihood(d, m, sP, rows, columns);

			if (u<sll+llP-ll){ // If accepted
//			if (u<llP-ll){ // If accepted
//				printf("Sigma: llP-ll = %g - %g\n",llP,ll);
//				printf("Accepted!\n");
				ll=llP;
				copyArray(sP,columns,s);
//				for (int n=0; n<np+2; n++){
//					printf("%u,",bndPntsP[n]);
//				}
//				printf("\n%u\n",np);
//				printf("sigma: ");
//				for (int n=0; n<columns; n++){
//					printf("%g,",s[n]);
//				}
//				printf("\n");
			}
		} 


		// Add a changepoint
		else if (r < MOVE+SIGMA+BIRTH && np<npmax){
			copyArrayUint(bndPntsC,K+2,bndPntsP);

			// Pick which changepoint to add right of
			pick=pcg32_random_r(&rng)/(RAND_MAX_U32/(np+1));
			bndPntsP[np+2] = pcg32_random_r(&rng) / ( RAND_MAX_U32 / ((bndPntsP[pick+1]-1)-(bndPntsP[pick]+1)+1) ) + (bndPntsP[pick]+1);
			npP=unique_uints(bndPntsP,np+3)-2;

			// Update the model
			update_model(bndPntsP, npP, rows, columns, &rng, d, s, m);

			// Calculate log likelihood for proposal
			lqz=lqxz(d, m, &bndPntsP[pick], rows, columns);
			llP=log_likelihood(d, m, s, rows, columns);
			if (u<-lqz+llP-ll){
//			if (u<llP-ll){
//				printf("Birth: llP-ll = %g - %g\n",llP,ll);
//				printf("Accepted!\n");
				ll=llP;
				np=npP;
				copyArrayUint(bndPntsP,K+2,bndPntsC);
//				printf("Birth, lqz = %g\n",-lqz);
				for (int n=0; n<np; n++){
					printf("%u,",bndPntsP[n+1]);
				}
				printf("\b\n");
//				printf("\n%u\n",np);
			}
		} 
		

		// Delete a changepoint
		else if (r < MOVE+SIGMA+BIRTH+DEATH && np>npmin){
			copyArrayUint(bndPntsC,K+2,bndPntsP);

			// Pick which changepoint to delete
			pick=pcg32_random_r(&rng)/(RAND_MAX_U32/np)+1;
			bndPntsP[pick]=bndPntsP[pick+1];
			npP=unique_uints(bndPntsP,np+2)-2;

			// Update the model
			update_model(bndPntsP, np, rows, columns, &rng, d, s, m);

			// Calculate log likelihood for proposal
			llP=log_likelihood(d, m, s, rows, columns);
			lqz=lqxz(d, m, &bndPntsC[pick], rows, columns);
			if (u<lqz+llP-ll){
//			if (u<llP-ll){
//				printf("Death: llP-ll = %g - %g\n",llP,ll);
//				printf("Accepted!\n");
				ll=llP;
				np=npP;
				copyArrayUint(bndPntsP,K+2,bndPntsC);
//				printf("Death, lqz = %g\n",lqz);
				for (int n=0; n<np; n++){
					printf("%u,",bndPntsP[n+1]);
				}
				printf("\b\n");
//				printf("\n%u\n",np);
			}
		} 
		

		// Update the model
		else {
			// Update the model
			update_model(bndPntsP, np, rows, columns, &rng, d, s, m);

			// Calculate log likelihood for proposal
			llP=log_likelihood(d, m, s, rows, columns);
			if (u<llP-ll){
//				printf("Update: llP-ll = %g - %g\n",llP,ll);
//				printf("Accepted!\n");
				ll=llP;
//				printf("Death, lqz = %g\n",lqz);
//				for (int n=0; n<np+2; n++){
//					printf("%u,",bndPntsP[n]);
//				}
//				printf("\n%u\n",np);

			}
		}
//		printf("%u\n",i);


	}

	return 0;
}


