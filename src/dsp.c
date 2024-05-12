#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>

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

// Define os símbolos da constelação 4QAM
const float complex QAM4_symbols[4] = {
    1.0 + 1.0*I, // 01
    1.0 - 1.0*I, // 11
   -1.0 + 1.0*I, // 00
   -1.0 - 1.0*I  // 10
};

/*
    Mapeia uma sequência de bits para símbolos da constelação 4-QAM.
Input: 
    bits (int): sequência de bits pseudo-aleatória.
    length (int): comprimento do array.
Output: 
    symbols (float complex): apontador para um vetor de símbolos mapeados.
*/

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
    Reduz a taxa de amostragem do sinal.
Input: 
    signal (float complex): sinal de entrada para decimação.
    length (int): comprimento do sinal.
    factor (int): fator de downsampling.
Output: 
    signalDown (float complex): apontador para um sinal decimado.
*/

float complex *downsample(float complex* signal, int length, int factor) {
    int newLength = length / factor;
    float complex *signalDown = (float complex*)malloc(newLength * sizeof(float complex));
    
    for(int i = 0; i < newLength; i++) {
        signalDown[i] = signal[i * factor];
    }
    
    return signalDown;
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

/*
    Executa a filtragem FIR e compensa o atraso do filtro.
Input: 
    x (float complex): sinal a ser convoluído.
    h (float complex): coeficientes do filtro.
    length_x (int): comprimento do sinal x.
    length_h (int): comprimento do sinal h
Output: 
    y (float complex): sinal filtrado.
*/

float complex *firFilter(float complex* x, float complex* h, int length_x, int length_h) {
    float complex *y = (float complex*)malloc(length_x * sizeof(float complex));
    // Calcula o índice de deslocamento para a convolução "same"
    int shift = (length_h - 1) / 2;

    for(int n = 0; n < length_x; n++) {
        // Inicializa o filtro com zeros
        y[n] = 0.0 + 0.0*I;
        for(int k = 0; k < length_x; k++) {
            if((n - k + shift) >= 0 && (n - k + shift) < length_h) {
                y[n] += x[k] * h[n - k + shift];
            }
        }
    }
    return y;
}


float complex* matchedFilter(float complex* x, float complex* h, int length_x, int length_h){
    float complex *sigRx = (float complex*)malloc(length_x * sizeof(float complex));
    sigRx = firFilter(x, h, length_x, length_h);
    sigRx = pnorm(sigRx, length_x);
    return sigRx;
}

/*
    Função responsável pela simulação da geração de bits,
    modulação, normalização e upsampling.
Input: 
    SpS (int): amostras por símbolo.
    Nbits (int): comprimento do bitstream
Output: 
    symbolsUp (float complex): apontador para a sequência de bits mapeados após upsampling.
*/

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

/*
    Retorna os índices dos valores mínimos.
Input: 
    array (float): vetor |r-sm|**2.
    length (int): comprimento do array
Output: 
    minIndex (int): índice do símbolo da constelação com a menor métrica de distância.
*/

int argmin(float *array, int length) {
    if (length <= 0) {
        // Se o comprimento for inválido, retorne -1
        return -1;
    }

    int minIndex = 0;
    float minValue = array[0];

    for (int i = 1; i < length; i++) {
        if (array[i] < minValue) {
            // Se o elemento atual for menor que o valor mínimo atual,
            // atualize o índice mínimo e o valor mínimo
            minIndex = i;
            minValue = array[i];
        }
    }

    return minIndex;
}

/*
    Execute a detecção de símbolos usando a regra ML (Máxima Verossimilhança)
    - Maximum Likelihood rule
Input: 
    r (float complex): sinal recebido.
    constSymb (float complex): os símbolos da constelação.
    length_constSymb (int) comprimento da matriz de símbolos.
    length_r (int) comprimento do sinal recebido.
Output: 
    decided (float complex): os símbolos detectados.
*/

float complex* MLdetector(float complex* r, float complex* constSymb, int length_constSymb, int length_r){

    float complex *decided = (float complex*)malloc(length_r * sizeof(float complex));
    // aloca a matriz de indexação para os símbolos
    int *indDec = (int*)malloc(length_r * sizeof(int));
    // matriz de valores |r-sm|**2, para m = 1,2,...,M
    float *distMetric = (float*)malloc(length_constSymb * sizeof(float));

    for(int ri = 0; ri < length_r; ri++){
        for(int index = 0; index < length_constSymb; index++){
            // calcula |r-sm|**2, para m = 1,2,...,M 
            distMetric[index] = pow(cabsf(r[ri] - constSymb[index]), 2);
        }
        // encontre o símbolo da constelação com a menor métrica de distância.
        indDec[ri] = argmin(distMetric, length_constSymb);
        // tome a decisão em favor do símbolo com a menor métrica
        decided[ri] = constSymb[indDec[ri]];
    }
    return decided;
}