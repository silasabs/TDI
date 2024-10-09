# Projeto final da disciplina de TDI

Abaixo encontra-se o fluxo da simulação baseada em código C para uma transmissão digital utilizando 4-QAM. A demodulação e realizada utilizando o critério de decisão ótima para o caso em que os símbolos enviados pelo transmissor são equiprováveis (*critério de máxima verossimilhança*).

<p align="center">
<img src=https://i.postimg.cc/Vv7zyh2p/Screenshot-from-2024-10-08-22-03-37.png>
</p>

## Recursos Implementados

- Transformação de Box-Muller
- Ruído Gaussiano 
- Mapeador Digital 4-QAM
- Upsample / Downsample
- Filtragem FIR / Correspondente
- Detector de Máxima Verossimilhança

O projeto utiliza o ctypes para chamar funções escritas em linguagens de baixo nível, como o C, diretamente a partir do Python. Para obter uma DLL a partir do GCC:

Linux Distros
```
$ gcc -shared -o dsp.so dsp.c
```

Microsoft Windows
```
$ gcc -shared -o dsp.dll dsp.c
```