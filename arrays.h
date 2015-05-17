#include <stdlib.h>
#include <math.h>

// Make a 1d array of points spaced by Step between Lower and Upper
double* array(double lower, double step, double upper){
	if (step<=0){
		perror("Step size too small");
	} 
	else {
		int i=0, numsteps=(upper-lower)/step+1;
		double *someArray=malloc(numsteps*sizeof(double));
		while (i<numsteps) {
			someArray[i]=lower+i*step;
			i++;
		}
		return someArray;
	}
	return NULL;
}

int* arrayInt(int lower, int step, int upper){
	if (step<=0){
		perror("Step size too small");
	} 
	else {
		int i=0, numsteps=(upper-lower)/step+1;
		int  *someArray=malloc(numsteps*sizeof(int));
		while (i<numsteps) {
			someArray[i]=lower+i*step;
			i++;
		}
		return someArray;
	}
	return NULL;
}

/* Make a 1-d array of n points linearly spaced between Lower and Upper */
double* linspace(double lower, double upper, double n){
	int i=0;
	double *someArray=malloc(n*sizeof(double));
	while (i<n){
		someArray[i]=(upper-lower)*i/(n-1)+lower;
		i++;
	}
	return someArray;
}

int* linspaceInt(int lower, int upper, int n){
	int i=0, *someArray=malloc(n*sizeof(int));
	while (i<n){
		someArray[i]=(upper-lower)*i/(n-1)+lower;
		i++;
	}
	return someArray;
}

/* Find the minimum of a 1-d array */
double min(double a, double b){
	if (b<a) return b;
	return a;
}

int minInt(int a, int b){
	if (b<a) return b;
	return a;
}

double minArray(double someArray[],int length){
	int i;
	double currentMin=someArray[0];
	for (i=1; i<length; i++){
		if (someArray[i]<currentMin){currentMin=someArray[i];}
	}
	return currentMin;
}

int minArrayInt(int someArray[],int length){
	int i, currentMin=someArray[0];
	for (i=1; i<length; i++){
		if (someArray[i]<currentMin){currentMin=someArray[i];}
	}
	return currentMin;
}

/* Find the maximum of a 1-d array */
double max(double a, double b){
	if (b>a) return b;
	return a;
}

int maxInt(int a, int b){
	if (b>a) return b;
	return a;
}

double maxArray(double someArray[], int length){
	int i;
	double currentMax=someArray[0];
	for (i=1; i<length; i++){
		if (someArray[i]>currentMax){currentMax=someArray[i];}
	}
	return currentMax;
}

int maxArrayInt(int someArray[], int length){
	int i, currentMax=someArray[0];
	for (i=1; i<length; i++){
		if (someArray[i]>currentMax){currentMax=someArray[i];}
	}
	return currentMax;
}



/* Malloc a 2D array of doubles of dimensions Rows x Columns */
double **mallocDoubleArray(int rows, int columns) {
	int i; double **someArray;
	someArray = (double **) malloc(rows*sizeof(double *));
	for (i = 0; i < rows; i++){
		someArray[i] = (double *) malloc(columns*sizeof(double));
	}
	return someArray;
} 

/* Free a 2D array of doubles of length rows */
void freeDoubleArray(double **someArray, int rows) {
	int i;
	for (i=0; i<rows; i++) {free(someArray[i]);}
	free(someArray);
}


/* Malloc a 2D array of ints of dimensions Rows x Columns */
int **mallocIntArray(int rows, int columns) {
	int i; int **someArray;
	someArray = (int **) malloc(rows*sizeof(int *));
	for (i = 0; i < rows; i++){
		someArray[i] = (int *) malloc(columns*sizeof(int));
	}
	return someArray;
} 

/* Free a 2D array of ints of length rows */
void freeIntArray(int **someArray, int rows) {
	int i;
	for (i=0; i<rows; i++) {free(someArray[i]);}
	free(someArray);
}




/* Parse a CSV file, return pointer to a 2d array of doubles */
double **csvparse(const char filePath[], const char delim, int * const restrict rows, int * const restrict columns){

	// File to open
	FILE *fp;
	fp=fopen(filePath,"r");

	if (fp==NULL){
		printf("WARNING: file does not exist!\n");
		return NULL;

	} else {

		char c;
		int numColumns=0, maxColumns=0, numChars=0, maxChars=0, numRows=0;

		// Determine maximum number of characters per row, delimiters per row, and rows
		for(c=getc(fp); c != EOF; c=getc(fp)){
			numChars++;
			if (c==delim){numColumns++;}
			if (c=='\n'){
				numRows++;
				//if there is a trailing delimiter, don't add an extra column for it
				fseek(fp,-2,SEEK_CUR);
				if(getc(fp)!=delim){numColumns++;}
				fseek(fp,+1,SEEK_CUR);	
				// See if we have a new maximum, and reset the counters
				if (numChars>maxChars){maxChars=numChars;}
				if (numColumns>maxColumns){maxColumns=numColumns;}
				numChars=0;
				numColumns=0;
			}
		} 
		// If the last line isn't blank, add one more to the row counter
		fseek(fp,-1,SEEK_CUR);
		if(getc(fp)!='\n'){numRows++;}


		// Malloc space for the imported array
		double ** const restrict importedMatrix=mallocDoubleArray(numRows,maxColumns);
		*rows=numRows;
		*columns=maxColumns;				

		// For each line,
		int i=0, j=0, k, field[maxColumns+2];
		char str[maxChars+2];
		char *endp;
		rewind(fp);
		while(fgets(str,maxChars+2,fp)!=NULL){

			// identify the delimited fields,
			field[0]=0;
			for(k=1;k<maxColumns+2;k++){field[k]=maxChars+1;}
			for(k=0, numColumns=0; str[k]!='\0'; k++){
				if (str[k]==delim){
					str[k]='\0';
					field[numColumns+1]=k+1;
					numColumns++;
				} else if(str[k]=='\n'){str[k]='\0';}
			}

			// and perform operations on each field
			for (j=0;j<maxColumns;j++){
				importedMatrix[i][j]=strtod(&str[field[j]],&endp);
				if(endp==&str[field[j]]){importedMatrix[i][j]=NAN;}
			}
			i++;
		}
		fclose(fp);

		printf("Maximum number of characters: %d\n", maxChars);
		printf("Maximum number of delimiters: %d\n", maxColumns);
		printf("Number of rows: %d\n", numRows);

		return importedMatrix;
	}

}


/* Parse a CSV file, return pointer to a 2d array of doubles */
double *csvparseflat(const char filePath[], const char delim, int * const restrict rows, int * const restrict columns){

	// File to open
	FILE *fp;
	fp=fopen(filePath,"r");

	if (fp==NULL){
		printf("WARNING: file does not exist!\n");
		return NULL;

	} else {

		char c, cl;
		int numColumns=0, maxColumns=0, numChars=0, maxChars=0, numRows=0;

		// Determine maximum number of characters per row, delimiters per row, and rows
		for(c=getc(fp); c != EOF; c=getc(fp)){
			numChars++;
			if (c==delim){numColumns++;}
			if (c=='\n'){
				numRows++;
				numColumns++;
				fseek(fp,-2,SEEK_CUR);
				cl=getc(fp);
				//if there is a trailing delimiter, don't add an extra column for it
				if (cl==delim){
					numColumns--;
				//if there is an empty line, don't add an extra row for it
				} else if (cl=='\n'){
					numRows--;
				}
				fseek(fp,+1,SEEK_CUR);	
				// See if we have a new maximum, and reset the counters
				if (numChars>maxChars){maxChars=numChars;}
				if (numColumns>maxColumns){maxColumns=numColumns;}
				numChars=0;
				numColumns=0;
			}
		} 
		// If the last line isn't blank, add one more to the row counter
		fseek(fp,-1,SEEK_CUR);
		if(getc(fp)!='\n'){numRows++;}


		// Malloc space for the imported array
		double * const restrict importedMatrix=malloc(sizeof(double)*maxColumns*numRows);
		*rows=numRows;
		*columns=maxColumns;				

		// Read the array
		int i=0, j=0, k, field[maxColumns+2];
		char str[maxChars+2];
		char *endp;
		rewind(fp);
		// For each line,
		while(fgets(str,maxChars+2,fp)!=NULL){
			// if line is not empty,
			if (str[0]!='\n'){
				// identify the delimited fields,
				field[0]=0;
				for(k=1;k<maxColumns+2;k++){field[k]=maxChars+1;}
				for(k=0, numColumns=0; str[k]!='\0'; k++){
					if (str[k]==delim){
						str[k]='\0';
						field[numColumns+1]=k+1;
						numColumns++;
					} else if(str[k]=='\n'){str[k]='\0';}
				}

				// perform operations on each field
				for (j=0;j<maxColumns;j++){
					importedMatrix[j*numRows + i]=strtod(&str[field[j]],&endp);
					if(endp==&str[field[j]]){importedMatrix[j*numRows + i]=NAN;}
				}
				// and increment the row counter
				i++;
			}
		}
		fclose(fp);

		printf("Maximum number of characters: %d\n", maxChars);
		printf("Maximum number of delimiters: %d\n", maxColumns);
		printf("Number of rows: %d\n", numRows);

		return importedMatrix;
	}

}

int csvwriteflat(const double* data, const char filePath[], const char* mode,  const char delim, int rows, int  columns){
	// File to open
	FILE *fp;
	fp=fopen(filePath, mode);

	//Print the array
	for (int i=0;i<rows;i++){
		for (int j=0;j<columns;j++){
			fprintf(fp, "%g%c", data[j*rows + i], delim);
		}
		// Delete the trailing delimiter and print a newline at the end of the row
		fseek(fp,-1,SEEK_CUR);
		fprintf(fp,"\n");
	}
	fclose(fp);
	return 0;
}


// Compute mean and variance with Knuth's algorithm
void Knuth_var(double* x, int n, double* mean, double* var){
	double delta, mu=0, m2=0;

	for (int i=0; i<n; i++){
		delta = x[i] - mu;
		mu += delta / (i+1);
		m2 += delta*(x[i] - mu);
	}

	if (n>1){
		*mean = mu;
		*var = m2/(n-1);
	} else {
		*mean = x[0];
		*var = 0;
	}
}

// Compute mean and variance with Knuth's algorithm, ignoring NaN data
void nanvar(double* x, int n, double* mean, double* var){
	double delta, mu=0, m2=0;
	int exists=0;

	for (int i=0; i<n; i++){
		if (!isnan(x[i])){
			exists++;
			delta = x[i] - mu;
			mu += delta / exists;
			m2 += delta*(x[i] - mu);
		}
	}

	if (exists>1){
		*mean = mu;
		*var = m2/(exists-1);
	} else if (exists==1 && n>1){
		*mean = mu;
		*var = 0;
	} else {
		*mean = NAN;
		*var = NAN;
	}
}

// Compute mean and standard deviation with Knuth's algorithm, ignoring NaN data
void nanstd(double* x, int n, double* mean, double* std){
	double delta, mu=0, m2=0;
	int exists=0;

	for (int i=0; i<n; i++){
		if (!isnan(x[i])){
			exists++;
			delta = x[i] - mu;
			mu += delta / exists;
			m2 += delta*(x[i] - mu);
		}
	}

	if (exists>1){
		*mean = mu;
		*std = sqrt(m2/(exists-1));
	} else if (exists==1 && n>1){
		*mean = mu;
		*std = 0;
	} else {
		*mean = NAN;
		*std = NAN;
	}
}


// Compute sum of a double array, excluding NaNs
double nansum(double* x, int n){
	double S=0.0;
	int exists=0;
	for (int i=0; i<n; i++){
		if (!isnan(x[i])){
			S+=x[i];
			exists=1;
		}
	}
	if (exists){
		return S;
	} else {
		return NAN;
	}
}

// Compute sum of the squares of a double array, excluding NaNs
double nansumSq(double* x, int n){
	double S=0.0;
	int exists=0;
	for (int i=0; i<n; i++){
		if (!isnan(x[i])){
			S+=x[i]*x[i];
			exists=1;
		}
	}
	if (exists){
		return S;
	} else {
		return NAN;
	}

}


// Compute mean of a double array, excluding NaNs
double nanmean(double* x, int n){
	double S=0.0;
	int exists=0;
	for (int i=0; i<n; i++){
		if (!isnan(x[i])){
			exists++;
			S+=x[i];
		}
	}
	return S/exists;
}

// Normalize a double array, excluding NaNs
void normalize(double* x, int n){
	double mu=nanmean(x,n);
	for (int i=0; i<n; i++){
		x[i]-=mu;
	}
}


// From wikipedia (public domain)
/* Comparison function. Receives two generic (void) pointers. */
int compare_ints(const void *p, const void *q){
	int x = *(const int *)p;
	int y = *(const int *)q;
	/* to avoid undefined behaviour through signed integer overflow,
	 *         avoid: return x - y; */
	int ret;
	if (x == y){
		ret = 0;
	} else {
		ret = (x < y) ? -1 : 1;
	}
	return ret;
}
 
/* Sort an array of n integers, pointed to by a. */
void sort_ints(int *a, size_t n){
	qsort(a, n, sizeof(int), compare_ints);
}

