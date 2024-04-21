#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>

/*
    Normaliza a potência média de cada componente de nums.
Input: 
    nums (float complex): fasores da constelação.
    size (int): comprimento do sinal
Output: 
    normalized_nums (float complex): Sinal com cada componente normalizado em potência.
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
    M (int): Ordem do esquema de modulação.
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

int main() {

    int tamanho = 3; // Tamanho da matriz

    float complex* numeros = malloc(tamanho * sizeof(float complex));

    // Preenche 'numeros' com os valores especificados

    numeros[0] = 1.0 + 2.0 * I;

    numeros[1] = 2.0 - 3.0 * I;

    numeros[2] = -1.0 - 1.0 * I;

    // Chama a função pnorm para normalizar a matriz

    float complex* resultado = pnorm(numeros, tamanho);

    for(int i = 0; i < tamanho; i++){

        printf("%.2f + %.2fi\n", creal(resultado[i]), cimag(resultado[i]));

    }

    return 0;

}