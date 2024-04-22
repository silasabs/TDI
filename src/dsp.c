#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>

/*
    Aumente a resolução de um sinal inserindo zeros entre as amostras.
Input: 
    signal (float complex): sinal de entrada para aumentar a resolução.
    length (int): comprimento do sinal.
    factor (int): fator de upsampling. O sinal será aumentado inserindo
                  "factor - 1" zeros entre cada amostra original.
Output: 
    signalUp (float complex): apontador para um sinal ampliado com zeros inseridos entre as amostras.
*/

float complex *upsample(float complex* signal, int length, int factor) {
    float complex *signalUp = (float complex*)malloc(length * factor * sizeof(float complex));
    for(int i = 0; i < length; i++) {
        signalUp[i*factor] = signal[i];
        for(int j = 1; j < factor; j++) {
            signalUp[i*factor+j] = 0;
        }
    }
    return signalUp;
}

/*
    Normaliza a potência média de cada componente de nums.
Input: 
    nums (float complex): fasores da constelação.
    size (int): comprimento do sinal
Output: 
    normalized_nums (float complex): sinal com cada componente normalizado em potência.
*/

float complex *pnorm(float complex* nums, int size) {
    // parameters
    float complex sum = 0;
    // normalized constellation matrix
    float complex *normalized_nums = (float complex*)malloc(size * sizeof(float complex));
    // calculate the average of the phasors
    for(int i = 0; i < size; i++) {
        sum += nums[i]*conj(nums[i]);
    }
    // normalizes the QAM constellation
    for(int i = 0; i < size; i++) {
        normalized_nums[i] = nums[i]/sqrt(creal(sum/size));
    }
    return normalized_nums;
}

/*
    Implementa a geração de uma sequência de bits pseudo-aleatórios 
    que chegam ao transmissor.
Input: 
    M (int): ordem do esquema de modulação.
    Nbits (int): comprimento do bitstream
Output: 
    bits (int): sequência de bits pseudo-aleatórios.
*/

int *getRandomBits(int Nbits) {

    int size = Nbits;
    int *bits = (int *)malloc(size * sizeof(int));
    // initialize the random number generator
    srand(time(NULL));

    for(int indBits = 0; indBits < size; indBits++){
        bits[indBits] = rand() % 2;
    }
    return bits;
}