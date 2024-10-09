# Projeto final da disciplina de TDI

Abaixo encontra-se o fluxo da simulação baseada em código C para uma transmissão digital utilizando 4-QAM. 

A demodulação e realizada utilizando o critério de decisão ótima para o caso em que os símbolos enviados pelo transmissor são equiprováveis (*critério de máxima verossimilhança*).

<br>

<p align="center">
<img src=https://i.postimg.cc/Vv7zyh2p/Screenshot-from-2024-10-08-22-03-37.png>
</p>

## Recursos Implementados

- Transformação de Box-Muller
- Ruído Gaussiano 
- Mapeador Digital 4-QAM
- Upsample
- Downsample
- Filtragem FIR
- Filtragem Correspondente
- Detector de Máxima Verossimilhança

Para obter uma DLL compartilhada a partir do GCC faça:
```
$ gcc -shared -o dsp.so dsp.c
```