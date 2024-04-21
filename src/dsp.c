#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*
  Implementa a geração de uma sequência de bits pseudo-aleatórios 
  que chegam ao transmissor.
Input: 
    M (int): Ordem do esquema de modulação.
    Nbits (int): comprimento do bitstream
Output: 
*/

int *getRandomBits(int Nbits, int M) {

    int size = Nbits * log2(M);
    int *bits = (int *)malloc(size * sizeof(int));
    // initialize the random number generator
    srand(time(NULL));

    for(int indBits = 0; indBits < size; indBits++){
        bits[indBits] = rand() % 2;
    }
    return bits;
}