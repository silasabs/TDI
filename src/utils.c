// Autor: Silas João Bezerra Soares.
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
    Autor: Moisés Oliveira
    Retorna o valor da função sinc com argumento x.
Input: 
    x (float): Argumento.
Output: 
     Sinc (float): Valor da função sinc cujo argumento é x.
*/

float sinc(float x){
    return sin(x * pi) / (x * pi);
}

/*  Autor: Silas João Bezerra Soares.
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

float complex abs_max(float complex *vetor, int N){
    float abs_max = fabs(vetor[0]);
    for (int i = 1; i < N; i++){
        if (fabs(vetor[i]) > abs_max)
            abs_max = fabs(vetor[i]);
    }
    return abs_max;
}