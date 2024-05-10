#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>

// Define os símbolos da constelação 4QAM
const float complex QAM4_symbols[4] = {
    1.0 + 1.0*I, // 01
    1.0 - 1.0*I, // 11
   -1.0 + 1.0*I, // 00
   -1.0 - 1.0*I  // 10
};

float complex *qam4Mapper(int *bits, int length) {
    float complex *symbols = (float complex*)malloc(length * sizeof(float complex));
    int i, symbol_index = 0;
    for (i = 0; i < length; i += 2) {
        
        int bit1 = bits[i];
        int bit2 = bits[i + 1];

        if (bit1 == 0 && bit2 == 0) {
            symbols[symbol_index++] = QAM4_symbols[2]; // 00
        } else if (bit1 == 0 && bit2 == 1) {
            symbols[symbol_index++] = QAM4_symbols[0]; // 01
        } else if (bit1 == 1 && bit2 == 0) {
            symbols[symbol_index++] = QAM4_symbols[3]; // 10
        } else if (bit1 == 1 && bit2 == 1) {
            symbols[symbol_index++] = QAM4_symbols[1]; // 11
        }
    }
    return symbols;
}

float **complex2float(float complex *arr, int length){
    // Aloca espaço para a matriz de partes reais e imaginárias
    float **separated_nums = (float **)malloc(2*length * sizeof(float));

    separated_nums[0] = (float *)malloc(length * sizeof(float)); // Parte real
    separated_nums[1] = (float *)malloc(length * sizeof(float)); // Parte imaginária

    // Preenche a matriz com as partes reais e imaginárias
    for(int i = 0; i < length; i++) {
        separated_nums[0][i] = crealf(arr[i]); 
        separated_nums[1][i] = cimagf(arr[i]); 
    }
    return separated_nums;
}

/*
    Leva uma sequência de bits aleatória para código gray.
Input: 
    bits (int): sequência de bits pseudo-aleatória.
    length (int): comprimento do array.
Output: 
    gray (int): apontador para uma sequência mapeada em código gray.
*/

int *grayMapping(int bits[], int length){
    int *gray = (int *)malloc(length * sizeof(int));
    // o primeiro bit e inalterado.
    gray[0] = bits[0];
    // obtém a sequência em código gray
    for(int i = 1; i < length; i++){
        gray[i] = bits[i-1] ^ bits[i];
    }
    return gray;
}

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
            signalUp[i*factor+j] = 0.0 + 0.0*I;
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
    float complex sum = 0;
    // matriz de constelação normalizada
    float complex *normalized_nums = (float complex*)malloc(size * sizeof(float complex));
    // calcular a média dos fasores
    for(int i = 0; i < size; i++) {
        sum += nums[i]*conj(nums[i]);
    }
    // normaliza a constelação
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
    int *bits = (int *)malloc(Nbits * sizeof(int));
    // initialize the random number generator
    srand(time(NULL));

    for(int indBits = 0; indBits < Nbits; indBits++){
        bits[indBits] = rand() % 2;
    }
    return bits;
}

float complex *mainUpSymbols(int Nbits, int SpS) {
    // gera os bits de forma pseudo-aleatória
    int *Bits = getRandomBits(Nbits);
    // mapeia os bits para os símbolos da constelação 4-QAM
    float complex *symbTx = qam4Mapper(Bits, Nbits);
    // normaliza os simbolos mapeados
    float complex *symbTxNorm = pnorm(symbTx, Nbits/2);
    // upsampling
    float complex *symbolsUp = upsample(symbTxNorm, Nbits/2, SpS);
    return symbolsUp;
}


// int main(){
//     // Vetor de bits de exemplo
//     int Nbits = 17*log2(4);
//     int SpS   = 16;
//     // Obtendo os bits aleatórios
//     int *Bits = getRandomBits(Nbits);
    
//     for (int i = 0; i < Nbits; i++) {
//         printf("%d", Bits[i]);
//     }
//     printf("\n");

//     // Mapeia os bits para os símbolos da constelação 4QAM
//     float complex *symbTx = qam4Mapper(Bits, Nbits);

//     // Imprime os símbolos mapeados
//     printf("Símbolos mapeados para a constelação 4QAM:\n");
//     for (int i = 0; i < Nbits/2; i++) {
//        printf("%.1f + %.1fi\n", crealf(symbTx[i]), cimagf(symbTx[i]));
//     }
    
//     // normaliza os simbolos mapeados para energia unitária
//     float complex *symbTxNorm = pnorm(symbTx, Nbits/2);

//     printf("Símbolos Normalizados:\n");
//     for (int i = 0; i < Nbits/2; i++) {
//        printf("%.4f + %.4fi\n", crealf(symbTxNorm[i]), cimagf(symbTxNorm[i]));
//     }

//     // Realiza o Upsampling
//     float complex *symbolsUp = upsample(symbTxNorm, Nbits/2, SpS);

//     printf("Após Upsampling:\n");
//     for (int i = 0; i < Nbits/2*SpS; i++) {
//        printf("%.4f + %.4fi\n", crealf(symbolsUp[i]), cimagf(symbolsUp[i]));
//     }

//     return 0;
// }

// testa pnorm e arrays de 2D
// int main() {
//     int size = 4; 

//     float complex *nums = (float complex *)malloc(size * sizeof(float complex));
    
//     // Inicializa o vetor de números complexos (apenas para teste)
//     for(int i = 0; i < size; i++) {
//         nums[i] = i + 1 + I * (i + 2);
//     }

//     // Normaliza os números complexos
//     float complex *normalized_nums = pnorm(nums, size);

//     // Converte os números complexos normalizados para um vetor de partes reais e imaginárias
//     float **float_nums = complex2float(normalized_nums, size);
    
//     // Imprime os resultados
//     printf("Parte Real Normalizada:\n");
//     for(int i = 0; i < size; i++) {
//         printf("%f ", float_nums[0][i]);
//     }
//     printf("\nParte Imaginária Normalizada:\n");
//     for(int i = 0; i < size; i++) {
//         printf("%f ", float_nums[1][i]);
//     }

//     // Libera a memória alocada
//     free(nums);
//     free(normalized_nums);
//     free(float_nums);
    
//     return 0;
// }

// testa codigo gray
// int main(){
//     int Nbits = 11;

//     // Obtendo os bits aleatórios
//     int *binary = getRandomBits(Nbits);
    
//     for (int i = 0; i < Nbits; i++) {
//         printf("%d", binary[i]);
//     }
//     printf("\n");
//     // Chamando a função para converter para código Gray
//     int *graybin = grayMapping(binary, Nbits);
//     for (int i = 0; i < Nbits; i++) {
//         printf("%d", graybin[i]);
//     }
//     printf("\n");
//     // Liberando a memória alocada para os bits
//     free(binary);

//     return 0;
// }