/* This program implements the Mersenne twister algorithm for generation of pseudorandom numbers. 
The program returns random integers in the range 0 to 2^32-1 (this holds even if a long int is
larger than 32 bits). Timing with gcc indicates that it is about twice as fast as the built in 
rand function. The original code was written by Michael Brundage and has been placed in the 
public domain.  */

#define MT_LEN 624
#include <stdlib.h>

void mt_init(unsigned long* mt_buffer, int* mt_index) {
    int i;
    for (i = 0; i < MT_LEN; i++){
        mt_buffer[i] = random();
    }
    *mt_index = 0;
}

#define MT_IA           397
#define MT_IB           (MT_LEN - MT_IA)
#define UPPER_MASK      0x80000000
#define LOWER_MASK      0x7FFFFFFF
#define MATRIX_A        0x9908B0DF
#define TWIST(mt_buffer,i,j)    ((mt_buffer)[i] & UPPER_MASK) | ((mt_buffer)[j] & LOWER_MASK)
#define MAGIC(s)        (((s)&1)*MATRIX_A)

unsigned long mt_rand_r(unsigned long* mt_buffer, int* mt_index) {
//    unsigned long * b = mt_buffer;
    int idx = *mt_index;
    unsigned long s;
    int i;
	
    if (idx == MT_LEN*sizeof(unsigned long))
    {
        idx = 0;
        i = 0;
        for (; i < MT_IB; i++) {
            s = TWIST(mt_buffer, i, i+1);
            mt_buffer[i] = mt_buffer[i + MT_IA] ^ (s >> 1) ^ MAGIC(s);
        }
        for (; i < MT_LEN-1; i++) {
            s = TWIST(mt_buffer, i, i+1);
            mt_buffer[i] = mt_buffer[i - MT_IB] ^ (s >> 1) ^ MAGIC(s);
        }
        
        s = TWIST(mt_buffer, MT_LEN-1, 0);
        mt_buffer[MT_LEN-1] = mt_buffer[MT_IA-1] ^ (s >> 1) ^ MAGIC(s);
    }
    *mt_index = idx + sizeof(unsigned long);
    return *(unsigned long *)((unsigned char *) mt_buffer + idx);
}
