#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <stdbool.h>

#define PI 3.14159265359

/**
 * @brief Converte matrizes complexas em matrizes de duas 
 *        dimensões de valor real e imaginário.
 * 
 * @param arr: matriz de números complexos.
 * @param length: comprimento da matriz.
 * @return float**: apontador duplo para uma matriz de representação real e imaginária. 
 *
 * Autor: Silas João Bezerra Soares.
 */

float **complex2float(float complex *arr, int length){
    // Aloca espaço para a matriz de partes reais e imaginárias
    float **separated_nums = (float **)malloc(2*length * sizeof(float));

    separated_nums[0] = (float *)malloc(length * sizeof(float)); // Parte real
    separated_nums[1] = (float *)malloc(length * sizeof(float)); // Parte imaginária

    // Preenche a matriz com as partes reais e imaginárias
    for (int i = 0; i < length; i++) {
        separated_nums[0][i] = crealf(arr[i]); 
        separated_nums[1][i] = cimagf(arr[i]); 
    }
    return separated_nums;
}

float complex abs_max(float complex *vetor, int N){
    float abs_max = fabs(vetor[0]);
    for (int i = 1; i < N; i++){
        if (fabs(vetor[i]) > abs_max)
            abs_max = fabs(vetor[i]);
    }
    return abs_max;
}

/* 
    Autor: Moisés Oliveira
    Retorna o valor da função sinc com argumento x.
Input: 
    x (float): Argumento.
Output: 
     Sinc (float): Valor da função sinc cujo argumento é x.
*/

float sinc(float x) {
    return sin(x * PI) / (x * PI);
}

/*
    Autor: Moisés Oliveira
    Gera um vetor cujos elementos são igualmente espaçados.
Input: 
    start (float): Primeiro elemento do vetor.
    stop (float): último elemento do vetor.
Output: 
    arr (array): Apontador para intervalos igualmente espaçados. 
*/

float* linspace(float start, float stop, int num) {
    float *arr = (float*)malloc(num * sizeof(float));
    float step = (stop - start) / (num - 1);
    for (int i = 0; i < num; i++) {
        arr[i] = start + i * step;
    }
    return arr;
}

/*
    Autor: Moisés Oliveira
    Retorna os coeficientes do filtro RC (Raised Cosine);
Input:
    t (float): Apontador para um vetor de tempo.
    length_t (int): comprimento do vetor de tempo.
    alpha (float): parâmetro roll-off do pulso cosseno levantado.
    Ts (float): Período de sinalização.
Output: 
    coeffs (float): Apontador para os coeficientes do filtro.
*/

float complex *rcFilterTaps(float *t,int length_t, float alpha, float Ts){
    float complex *coeffs = (float complex*)malloc(length_t * sizeof(float complex));
    for (int i=0; i < length_t; i++) {
        if (abs(t[i]) == Ts/(2*alpha)) {
            coeffs[i] = PI / (4*Ts) * sinc(1 / (2*alpha));
        }
        else {
            coeffs[i] = 1/Ts*sinc(t[i]/Ts)*cos(alpha*PI*t[i]/Ts)/(1-4*pow(alpha, 2)*pow(t[i],2)/pow(Ts,2));
        }
    }
    return coeffs;
}

/*
    Autor: Moisés Oliveira
    Gera os coeficientes normalizados do filtro formatador de pulso.
Input:
    SpS (int): Amostars por símbolo.
    length_t (int): comprimento do vetor de tempo.
    alpha (float): parâmetro roll-off do pulso cosseno levantado.
    Ts (float): Período de sinalização.
Output: 
    filter_coeffs (float): Apontador para os coeficientes do filtro.
*/

float complex *pulseShape(int SpS, int N, float alpha , float Ts){
    float *t;
    float fa = (1/Ts)*SpS;

    float complex *filter_coeffs = (float complex*)malloc(N * sizeof(float complex));
    
    t = linspace(-N/2*1/fa, N/2*1/fa, N);
    filter_coeffs = rcFilterTaps(t, N, alpha, Ts);
    
    float sum_square = 0;
    for (int i = 0;i < N;i++) {
        sum_square = pow(filter_coeffs[i], 2) + sum_square;
    }
    for (int i = 0; i < N; i++) {
        filter_coeffs[i] = filter_coeffs[i] / sqrt(sum_square);
    }

    float max = abs_max(filter_coeffs, N);
    
    for (int i = 0; i < N; i++) {
        filter_coeffs[i] = filter_coeffs[i] / max;
    }
    
    return filter_coeffs;
}

float normal_random(float mu, float stddev) {

    float u1 = (float)rand() / RAND_MAX;
    float u2 = (float)rand() / RAND_MAX;

    float z = sqrt(-2 * log(u1)) * cos(2 * PI * u2);  
    return z * stddev + mu;
}

float complex* noise(int samples, float mu, float stddev) {
    
    float complex *noise = (float complex*)malloc(samples * sizeof(float complex));

    srand(time(NULL));

    for (int i = 0; i < samples; i++) {
        
        float real_part = normal_random(mu, stddev);
        float imag_part = normal_random(mu, stddev);
        
        noise[i] = real_part + I * imag_part;
    }

    return noise;
}

/**
 * @brief Define os símbolos da constelação 4QAM
 * Autor: Silas João Bezerra Soares.
 */

const float complex QAM4_symbols[4] = {
    1.0 + 1.0*I, // 01
    1.0 - 1.0*I, // 11
   -1.0 + 1.0*I, // 00
   -1.0 - 1.0*I  // 10
};

/**
 * @brief Mapeia uma sequência de bits para símbolos da constelação 4-QAM.
 * 
 * @param bits: sequência de bits pseudo-aleatória.
 * @param length: comprimento da sequência de bits.
 * @return float*: apontador para um vetor de símbolos mapeados.
 *
 * Autor: Silas João Bezerra Soares.
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

/**
 * @brief Leva uma sequência de bits aleatória para código gray.
 * 
 * @param bits: sequência de bits pseudo-aleatória.
 * @param length: comprimento do array.
 * @return int*: apontador para uma sequência mapeada em código gray. 
 *
 * Autor: Silas João Bezerra Soares.
 */

int *grayMapping(int bits[], int length){
    int *gray = (int *)malloc(length * sizeof(int));
    // o primeiro bit e inalterado.
    gray[0] = bits[0];
    // obtém a sequência em código gray
    for (int i = 1; i < length; i++) {
        gray[i] = bits[i-1] ^ bits[i];
    }
    return gray;
}

/**
 * @brief Aumente a resolução de um sinal inserindo zeros entre as amostras.
 * 
 * @param signal: sinal de entrada. 
 * @param length: comprimento do sinal. 
 * @param factor: fator de upsampling. O sinal será aumentado inserindo
                  "factor - 1" zeros entre cada amostra original.
 * @return float*: apontador para um sinal ampliado com zeros inseridos entre as amostras.
 *
 * Autor: Silas João Bezerra Soares.
 */

float complex *upsample(float complex* signal, int length, int factor) {
    float complex *signalUp = (float complex*)malloc(length * factor * sizeof(float complex));
    for (int i = 0; i < length; i++) {
        signalUp[i*factor] = signal[i];
        for (int j = 1; j < factor; j++) {
            signalUp[i*factor+j] = 0.0 + 0.0*I;
        }
    }
    return signalUp;
}

/**
 * @brief Reduz a taxa de amostragem do sinal.
 * 
 * @param signal: sinal de entrada. 
 * @param length: comprimento do sinal.
 * @param factor: fator de downsampling.
 * @return float*: apontador para um sinal decimado.
 *
 * Autor: Silas João Bezerra Soares.
 */

float complex *downsample(float complex* signal, int length, int factor) {
    int newLength = length / factor;
    float complex *signalDown = (float complex*)malloc(newLength * sizeof(float complex));
    
    for (int i = 0; i < newLength; i++) {
        signalDown[i] = signal[i * factor];
    }
    
    return signalDown;
}

/**
 * @brief Normaliza a potência média de cada componente de nums.
 * 
 * @param nums: fasores da constelação 
 * @param size: comprimento do sinal 
 * @return float*: sinal com cada componente normalizado em potência.
 *
 * Autor: Silas João Bezerra Soares.
 */

float complex *pnorm(float complex* nums, int size) {
    float complex sum = 0;
    // matriz de constelação normalizada
    float complex *normalized_nums = (float complex*)malloc(size * sizeof(float complex));
    // calcular a média dos fasores
    for (int i = 0; i < size; i++) {
        sum += nums[i]*conj(nums[i]);
    }
    // normaliza a constelação
    for (int i = 0; i < size; i++) {
        normalized_nums[i] = nums[i]/sqrt(creal(sum/size));
    }
    return normalized_nums;
}


/**
 * @brief Implementa a geração de uma sequência de bits pseudo-aleatórios 
 *        que chegam ao transmissor.
 * 
 * @param Nbits: comprimento do bitstream
 * @return int*: sequência de bits pseudo-aleatórios.
 *
 * Autor: Silas João Bezerra Soares.
 */

int *getRandomBits(int Nbits) {
    int *bits = (int *)malloc(Nbits * sizeof(int));
    // initialize the random number generator
    srand(time(NULL));

    for (int indBits = 0; indBits < Nbits; indBits++) {
        bits[indBits] = rand() % 2;
    }
    return bits;
}

/**
 * @brief Executa a filtragem FIR e compensa o atraso do filtro.
 * 
 * @param x: sinal a ser convoluído.
 * @param h: coeficientes do filtro.
 * @param length_x: comprimento do sinal x.
 * @param length_h: comprimento do sinal h.
 * @return float*: sinal filtrado.
 *
 * Autor: Silas João Bezerra Soares.
 */
 
float complex *firFilter(float complex* x, float complex* h, int length_x, int length_h) {
    float complex *y = (float complex*)malloc(length_x * sizeof(float complex));
    // Calcula o índice de deslocamento para a convolução "same"
    int shift = (length_h - 1) / 2;

    for (int n = 0; n < length_x; n++) {
        // Inicializa o filtro com zeros
        y[n] = 0.0 + 0.0*I;
        for (int k = 0; k < length_x; k++) {
            if ((n - k + shift) >= 0 && (n - k + shift) < length_h) {
                y[n] += x[k] * h[n - k + shift];
            }
        }
    }
    return y;
}

/**
 * @brief Aplica a filtragem correspondente.
 * 
 * @param x: sinal a ser convoluído.
 * @param h: coeficientes do filtro.
 * @param length_x: comprimento do sinal x.
 * @param length_h: comprimento da resposta ao impulso h
 * @return float*: sinal filtrado e normalizado.
 *
 * Autor: Silas João Bezerra Soares.
 */

float complex* matchedFilter(float complex* x, float complex* h, int length_x, int length_h){
    float complex *sigRx = (float complex*)malloc(length_x * sizeof(float complex));
    sigRx = firFilter(x, h, length_x, length_h);
    sigRx = pnorm(sigRx, length_x);
    return sigRx;
}

/**
 * @brief Função responsável pela simulação da geração de bits,
 *        modulação, normalização e upsampling.
 * 
 * @param Nbits: comprimento do bitstream
 * @param SpS: amostras por símbolo.
 * @return float*: apontador para a sequência de bits mapeados após upsampling.
 *
 * Autor: Silas João Bezerra Soares.
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

// Autor: Silas João Bezerra Soares.
// float complex* mainRx(float complex* x, float complex* h, int SpS, int length_x, int length_h, float complex* constSymb){
//     float complex *sigRx = (float complex*)malloc(length_x * sizeof(float complex));
//     // filtro casado
//     float complex* sigRx = firFilter(x, h, length_x, length_h);
//     float complex* sigRx = pnorm(sigRx, length_x);
//     // demodulação
//     //float complex *r = downsample(sigRx, length_x, SpS);
//     //float complex *decided = MLdetector(r, constSymb, length_constSymb, 4, length_x);
//     return sigRx;
// }

/*  Autor: Silas João Bezerra Soares.
    Função responsável pela simulação da geração de bits,
    modulação, normalização e upsampling.
Input: 
    SpS (int): amostras por símbolo.
    Nbits (int): comprimento do bitstream
Output: 
    symbolsUp (float complex): apontador para a sequência de bits mapeados após upsampling.
*/

float complex *mainTx(int Nbits, int SpS, int N, float alpha, float Ts) {
    // gera os bits de forma pseudo-aleatória
    int *Bits = getRandomBits(Nbits);
    // mapeia os bits para os símbolos da constelação 4-QAM
    float complex *symbTx = qam4Mapper(Bits, Nbits);
    // normaliza os simbolos mapeados
    float complex *symbTxNorm = pnorm(symbTx, Nbits/2);
    // upsampling
    float complex *symbolsUp = upsample(symbTxNorm, Nbits/2, SpS);
    // gera um formato de pulso RC (cosseno levantado)
    float complex *pulse = pulseShape(SpS, N, alpha, Ts);
    // filtro formatador de pulso
    float complex *sigTx = firFilter(symbolsUp, pulse, Nbits/2*SpS, N);
    // retorna o sinal em banda base
    return sigTx;
}

/**
 * @brief Retorna os índices dos valores mínimos.
 * 
 * @param array vetor: |r-sm|**2.
 * @param length: comprimento do array
 * @return int: índice do símbolo da constelação com a menor métrica de distância.
 *
 * Autor: Silas João Bezerra Soares.
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

/**
 * @brief Execute a detecção de símbolos usando a regra ML (Máxima Verossimilhança)
 *        - Maximum Likelihood rule
 * 
 * @param r: sinal recebido.
 * @param constSymb: os símbolos da constelação.
 * @param length_constSymb: comprimento da matriz de símbolos.
 * @param length_r: comprimento do sinal recebido.
 * @return float*: os símbolos detectados.
 *
 * Autor: Silas João Bezerra Soares.
 */

float complex* MLdetector(float complex* r, float complex* constSymb, int length_constSymb, int length_r){

    float complex *decided = (float complex*)malloc(length_r * sizeof(float complex));
    // aloca a matriz de indexação para os símbolos
    int *indDec = (int*)malloc(length_r * sizeof(int));
    // matriz de valores |r-sm|**2, para m = 1,2,...,M
    float *distMetric = (float*)malloc(length_constSymb * sizeof(float));

    for (int ri = 0; ri < length_r; ri++) {
        for (int index = 0; index < length_constSymb; index++) {
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

// int main() {
//     int length = 10;  // Tamanho do sinal
//     float mu = 0.0;   // Média do ruído gaussiano
//     float stddev = 0.05;  // Desvio padrão do ruído gaussiano

//     // Criar um sinal exemplo (pode ser zero ou um sinal real)
//     float complex *signal = (float complex*)malloc(length * sizeof(float complex));
//     for (int i = 0; i < length; i++) {
//         signal[i] = i + I * i;  // Exemplo: sinal complexo crescente
//     }

//     printf("Sinal original:\n");
//     for (int i = 0; i < length; i++) {
//         printf("signal[%d] = %f + %fi\n", i, creal(signal[i]), cimag(signal[i]));
//     }

//     // Chamar a função noise() para adicionar ruído ao sinal
//     signal = noise(signal, length, mu, stddev, true);

//     printf("\nSinal com ruído adicionado:\n");
//     for (int i = 0; i < length; i++) {
//         printf("signal[%d] = %f + %fi\n", i, creal(signal[i]), cimag(signal[i]));
//     }

//     // Gerar ruído gaussiano sem adicionar ao sinal (apenas ruído)
//     float complex *generated_noise = noise(signal, length, mu, stddev, false);

//     printf("\nRuído gerado (sem adicionar ao sinal):\n");
//     for (int i = 0; i < length; i++) {
//         printf("noise[%d] = %f + %fi\n", i, creal(generated_noise[i]), cimag(generated_noise[i]));
//     }

//     // Liberar memória
//     free(signal);
//     free(generated_noise);

//     return 0;
// }