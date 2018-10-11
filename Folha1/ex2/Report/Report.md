---
title: Folha 1, Ex. 2 --- Estutura Electrónica e Modelação Computacional
author: Miguel Murça, 2015235874
lang: pt-PT
documentclass: article
header-includes: |
    \usepackage{supertabular, listings, physics, multicol}
    \lstset{
        basicstyle=\small,
        frame=single
    }
    \makeatletter
    \def\fps@figure{h}
    \makeatother
---

# Energias de Ionização

Recorrendo ao programa de Química Quântica PSI4, começou por se calcular,
recorrendo ao método de Hartree-Fock na base STO-3G, a energia total das
moléculas de H<sub>2</sub> e HeH<sup>+</sup>para distâncias de ligação
variáveis (\ref{code:vardist}).

Os pontos `(d, E)` obtidos encontram-se apresentados em \ref{data:varH}
e \ref{data:varHeH}. A representação gráfica dos pontos encontra-se
nas figuras \ref{fig:varH} e \ref{fig:varHeH}

![Gráfico d(E) dos valores obtidos para a energia total da molécula de Hidrogénio,
	de acordo com o método de Hartree-Fock, na base STO-3G. Foi calculada a energia
	para 500 pontos, uniformemente espaçados entre 0.1\AA e 1.5\AA.
	\label{fig:varH}](bits/figures/ex2-H.png)

![Gráfico d(E) dos valores obtidos para a energia total da molécula HeH<sup>+</sup>,
	de acordo com o método de Hartree-Fock, na base STO-3G. Foi calculada a energia
	para 300 pontos, uniformemente espaçados entre 0.1\AA e 1.5\AA.
	\label{fig:varHeH}](bits/figures/ex2-HeH.png)

Simultaneamente ao cálculo dos pontos apresentados, identificaram-se também os pontos
de energia total mínima:

$$ \text{Distância de Energia Mínima}_{\text{H}_2} 	 = d_{\text{H}_2}	= \left ( 0.7116 \, , \, -1.1175 \right ) $$

$$ \text{Distância de Energia Mínima}_{\text{HeH}^+} = d_{\text{HeH}^+} = \left ( 0.9288 \, , \, -2.8522 \right ) $$

Tomando as estas distâncias, determinaram-se, recorrendo de novo ao programa PSI4,
as propriedades das moléculas onde a distância de ligação é igual à distância de energia mínima.

Os excertos relevantes do *output* obtido encontram-se apresentados nas figuras
\ref{output:Hmin} e \ref{output:HeHmin}

```{=latex}
\begin{figure}[htb]
\begin{lstlisting}
...
set basis sto-3g
...
energy('scf', molecule=H)
...
    Orbital Energies [Eh]
    ---------------------

    Doubly Occupied:                                                      

       1Ag    -0.590501  

    Virtual:                                                              

       1B1u    0.701093  

    Final Occupation by Irrep:
             Ag   B1g   B2g   B3g    Au   B1u   B2u   B3u 
    DOCC [     1,    0,    0,    0,    0,    0,    0,    0 ]

	Energy converged.

	@DF-RHF Final Energy:    -1.11752981987458
...
\end{lstlisting}
\caption{Excertos relevantes do output produzido pelo PSI4 para H2}
\label{output:Hmin}
\end{figure}
```

```{=latex}
\begin{figure}[htb]
\begin{lstlisting}
...
set basis sto-3g
...
energy('scf', molecule=HeH)
...
    Orbital Energies [Eh]
    ---------------------

    Doubly Occupied:                                                      

       1A1    -1.524298  

    Virtual:                                                              

       2A1    -0.267229  

    Final Occupation by Irrep:
             A1    A2    B1    B2 
    DOCC [     1,    0,    0,    0 ]


	Energy converged.

	@DF-RHF Final Energy:    -2.85441846621746
...
\end{lstlisting}
\caption{Excertos relevantes do output produzido pelo PSI4 para HeH+}
\label{output:HeHmin}
\end{figure}
```

Nestes estão indicadas as energias HOMO para cada espécie, nomeadamente:

$$ \epsilon^{\text{H}_2}_{\text{HOMO}} 	 = -0.590501 \quad \text{(Hartree)} $$

$$ \epsilon^{\text{HeH}^+}_{\text{HOMO}} = -1.270649 \quad \text{(Hartree)} $$

Assim como a energia total de cada espécie:

$$ \text{E}^{\text{H}_2} = -1.117530 \quad \text{(Hartree)} $$

$$ \text{E}^{\text{HeH}^+} = -2.854418 \quad \text{(Hartree)} $$

De acordo com o teorema de Koopmans, a energia de ionização $I_p$ de uma espécie química
com energia da HOMO $\epsilon_\text{HOMO}$ é dada por

$$ I_p = - \epsilon_\text{HOMO} $$

Este teorema assume que a contribuição da relaxação da orbital da espécie ionizada para
a energia total é desprezável face a $\epsilon$, ou seja, assume que as orbitais de
$\ket{{}^N \Psi_0}$ são óptimas para o estado $\ket{{}^{N-1} \Psi_a}$.
Esta aproximação pode ser designada por aproximação de "orbitais congeladas".
Desprezada a diminuição da energia pela relaxação, esperar-se-ia que a energia
prevista para a ionização da espécie fosse inferior (em módulo) à energia experimentalmente
verificada, mas é necessário ter em conta que se $\epsilon_\text{HOMO}$ for
obtido pelo método de Hartree-Fock, o termo de troca é perdido, pelo que a energia
obtida no cálculo direto da espécie ionizada (onde a relaxação é contabilizada)
pode não ser inferior (em módulo) à energia total obtida para a espécie não-ionizada.

Para verificar estas observações, calcularam-se, de novo recorrendo ao PSI4 com o método
de (*Unrestricted*) Hartree-Fock na base STO-3G (**UHF/STO-3G**), as energias totais para $\text{H}_2^+$ e
$\text{HeH}^{2+}$ (mantendo as respetivas distâncias de ligação internuclear
$d_{\text{H}_2}$ e $d_{\text{HeH}^+}$). Os ficheiros de input descrevendo as moléculas para
uso no PSI4 encontram-se apresentados na secção \ref{code:ioninp}.

Os excertos relevantes do *output* obtido encontram-se apresentados nas figuras
\ref{output:Hion}, \ref{output:HeHion}.

```{=latex}
\begin{figure}[htb]
\begin{lstlisting}
    Alpha Occupied:                                                       

       1Ag    -1.270649  

    Alpha Virtual:                                                        

       1B1u    0.032566  

    Beta Occupied:                                                        

    

    Beta Virtual:                                                         

       1Ag    -0.590501     1B1u    0.212115  

    Final Occupation by Irrep:
             Ag   B1g   B2g   B3g    Au   B1u   B2u   B3u 
    DOCC [     0,    0,    0,    0,    0,    0,    0,    0 ]
    SOCC [     1,    0,    0,    0,    0,    0,    0,    0 ]

  Energy converged.

  @DF-UHF Final Energy:    -0.52702876420193
...
\end{lstlisting}
\caption{Excertos relevantes do output produzido pelo PSI4 para $H_2^+$}
\label{output:Hion}
\end{figure}
```

```{=latex}
\begin{figure}[htb]
\begin{lstlisting}
    Alpha Occupied:                                                       

       1A1    -2.497832  

    Alpha Virtual:                                                        

       2A1    -0.837902  

    Beta Occupied:                                                        

    

    Beta Virtual:                                                         

       1A1    -1.480189     2A1    -0.751747  

    Final Occupation by Irrep:
             A1    A2    B1    B2 
    DOCC [     0,    0,    0,    0 ]
    SOCC [     1,    0,    0,    0 ]

  Energy converged.

  @DF-UHF Final Energy:    -1.35666548775990
...
\end{lstlisting}
\caption{Excertos relevantes do output produzido pelo PSI4 para $HeH^{2+}$}
\label{output:HeHion}
\end{figure}
```

A energia total das espécies ionizada é então, de acordo com UHF/STO-3G:

$$ \text{E}_{\text{H}_2^+} = -0.527029 \quad \text{(Hartree)} $$

$$ \text{E}_{\text{HeH}^{2+}} = -1.356665 \quad \text{(Hartree)} $$

Encontra-se, na tabela \ref{tab:ionenergies}, comparada a energia de ionização
obtida pelo teorema de Koopmans ($I_p$) com a
energia de ionização prevista por aplicação direta de UHF/STO-3G ($I_d$), assim como 
com energias de ionização experimentalmente previstas (NIST, 2018) ($I_e$).

|            	| $H_2$     	| $HeH^+$   	| $H_2^+$ 	| $HeH^{2+}$ 	|
|------------	|-----------	|---------------|---------------|------------	|
| E          	| -1.117553 	| -2.854418 	| -0.527029	| -1.356665     |
| $\epsilon$ 	| -0.590501    	| -1.270649 	| ---     	| ---        	|
| $I_p$      	|  0.590501   	|  1.270649  	| ---     	| ---        	|
| $I_d$      	|  0.590524  	|  1.497753  	| ---     	| ---        	|
| $I_e$         |  0.566895     |  0.432356[^1]	| ---       	| ---           |

Table:  Comparação das energias HOMO ($\epsilon$), total (E), e energias
        de ionização pelo teorema de Koopmans ($I_p$) e por cálculo da energia
        da espécie ionizada ($I_d$). Todos os valores estão expressos em Hartree.
        Os valores $E$ e $\epsilon$ foram obtidos por cálculo com UHF/STO-3G.
        \label{tab:ionenergies}

[^1]: O sistema HeH<sup>+</sup> é demasiado instável para ser determinada a energia de ionização
	experimentalmente. Como tal, tomou-se como valor de comparação o resultado obtido
	pelo método MP3=FULL/6-31+G** (perturbações de 2ª ordem)

\pagebreak

Observa-se que para o $H_2$ as energias de ionização $I_p$ e $I_d$ concordam em
quatro casas decimais, pelo que a energia de relaxação electrónica pode, de facto,
ser desprezada em boa aproximação. Isto é esperado, uma vez que o número de electrões
na espécie original ($H_2$) é muito reduzido, pelo que a energia de relaxação é pequena.

No entanto, observa-se que tanto $I_p$ como $I_d$ apenas concordam com $I_e$ na primeira
casa decimal; está presente um erro sistemático que poderá advir da escolha da base.
Recorrendo a UHF/cc-pVTZ (cujo input para o PSI4 se encontra apresentado em \ref{code:ccpvtz}
com output \ref{output:ccpvtz}), observam-se os resultados apresentados na tabela \ref{tab:ccpvtz}.

```{=latex}
\begin{minipage}{\textwidth}
```

|            	| $H_2$     	| $HeH^+$   	| $H_2^+$   	| $HeH^{2+}$ 	|
|------------	|-----------	|-----------	|-----------	|------------	|
| E          	| -1,132599 	| -2.923252 	| -0.560512 	| -1.447928  	|
| $\epsilon$ 	| -0.603595 	| -1.528449 	| ---       	| ---        	|
| $I_p$      	| 0.603595  	| 1.528449  	| ---       	| ---        	|
| $I_d$      	| 0.572087  	| 1.475324  	| ---       	| ---        	|

Table:  Valores obtidos para energia total (E), energia HOMO ($\epsilon$),
        energia de ionização a partir do teorema de Koopmans ($I_p$) e
        energia de ionização por cálculo da energia total da espécie ionizada
        ($I_d$), recorrendo a **URH/cc-pVTZ**. \label{tab:ccpvtz}

```{=latex}
\end{minipage}
```

Verifica-se que em comparação com a tabela \ref{tab:ionenergies},
a energia de ionização $I_d$ obtida encontra-se mais próxima da energia
experimentalmente verificada.

Para $\text{HeH}^+$ a carga nuclear é mais elevada, pelo que a energia de relaxação
deve ser mais elevada. Observa-se, de facto, que nesse caso $I_p$ e $I_d$ já são
discordantes na primeira casa decimal. Por outro lado, verifica-se que o Método
de Hartree-Fock não é adequado para o cálculo das energias de $\text{HeH}^+$;
ambas as bases produzem resultados com 1 Hartree de diferença do resultado de
referência.

# Energia de Dissociação

Para verificar qual a energia total dos sistemas (de acordo com UHF/STO-3G)
no regime de distância internuclear tão grande que se podem considerar átomos
separados, calcularam-se as energias totais para distâncias exponencialmente
crescentes, observando o valor energético assimptótico. Os resultados obtidos
encontram-se apresentados nas figuras \ref{fig:Hassimp}, \ref{fig:HeHassimp}.

![Valores de energia total da molécula H<sub>2</sub> (de acordo com UHF/STO-3G)
    para distâncias internucleares exponencialmente crescentes.
    \label{fig:Hassimp}](bits/figures/H.png)

![Valores de energia total da molécula de HeH<sup>+</sup> (de acordo com UHF/STO-3G)
    para distâncias internucleares exponencialmente crescentes.
    \label{fig:HeHassimp}](bits/figures/HeH+.png)

Observa-se que para distâncias da ordem de $10^4$ e superior as moléculas apresentam
comportamento assimptótico.

Neste regime, as energias totais das moléculas são (de acordo com UHF/STO-3G),

$$ \text{E}_{\text{H}_2}    = -0.546008 \quad \text{(Hartree)} $$

$$ \text{E}_{\text{HeH}^+}  = -2.807913 \quad \text{(Hartree)} $$

Calculou-se a energia dos átomos/iões individuais em UHF/STO-3G, obtendo-se
os resultados apresentados na tabela \ref{tab:atoms}.

+---------+--------------------------+
| Espécie | Energia Total (Hartree)  |
+---------+:------------------------:+
| H       | -0.466582                |
+---------+--------------------------+
| H+      | -0.000000                |
+---------+--------------------------+
| H-      | -0.158852                |
+---------+--------------------------+
| He      | -2.807913                |
+---------+--------------------------+
| He+     | -1.931748                |
+---------+--------------------------+

Table:  Valores obtidos (UHF/STO-3G) para a energia total
        dos átomos e respetivos iões nas moléculas de H<sub>2</sub>
        e HeH<sup>+</sup>. \label{tab:atoms}

Considerando as possíveis dissociações de H<sub>2</sub> e respetivas
energias (calculadas pela soma das energias das espécies dissociadas):

\clearpage
\begin{minipage}{\textwidth}

$\textbf{H}_2$
\vspace{1em}

\quad $E_\text{Grande Distância} = -0.546008 \quad \text{(Hartree)}$
\vspace{1em}

(A) \qquad $H_2 \longrightarrow H + H$

$$\sum E = -0.933164 \quad \text{(Hartree)}$$

(B) \qquad $H_2 \longrightarrow H^+ + H^-$

$$\sum E = -0.158852 \quad \text{(Hartree)}$$

\end{minipage}

\begin{minipage}{\textwidth}

$\textbf{HeH}^+$
\vspace{1em}

\quad $E_\text{Grande Distância} = -2.807913 \quad \text{(Hartree)}$
\vspace{1em}

(C) \qquad $HeH^{+} \longrightarrow He + H^+$

$$\sum E = -2.807913 \quad \text{(Hartree)}$$

(D) \qquad $HeH^{+} \longrightarrow He^+ + H$

$$\sum E = -2.398330 \quad \text{(Hartree)}$$

\end{minipage}
\vspace{1mm}

Prevê-se então que o $\text{HeH}^+$ se dissocie nas espécies
He e $\text{H}^+$, uma vez que a soma da energia total das
duas espécies resulta na energia assimptótica para $\text{HeH}^+$
quando a ligação internuclear tende para infinito.

Por outro lado, nenhuma das energias de (A) ou (B) corresponde à
energia total do sistema a grandes distâncias de ligação.
No entanto, tomando a média dos dois valores,

$$ \frac{\text{E}_\text{(A)} + \text{E}_\text{(B)}}{2} = -0.546008 \quad \text{(Hartree)} $$

verifica-se que o valor obtido coincide exatamente com a energia
total do sistema a grande distância internuclear.

Este resultado evidência a correlação introduzida no sistema pelo
método de Hartree-Fock; sendo sempre contabilizado o campo médio
(e sendo forçado o emparelhamento de spin em RHF),
a solução é dada pela sobreposição dos dois sistemas individuais
possíveis, e não pela solução física no caso
de distâncias tais que os dois átomos são independentes.

A previsão não é física, pois prevê que um dos eletrões não distinga
energeticamente os dois núcleos, independentemente da distância,
e portanto haja uma probabilidade 50/50 de dissociação em (A) e (B).

# Referências

NIST Computational Chemistry Comparison and Benchmark Database,\
NIST Standard Reference Database Number 101\
Release 19, April 2018, Editor: Russell D. Johnson III\
http://cccbdb.nist.gov/  

\topskip0pt
\vspace*{\fill}

----------------------------------------------------------------------------------------

\vspace*{\fill}
\appendix
\clearpage

# *Raw Data*

## H<sub>2</sub> 
\label{data:varH}

```{=latex}
\begin{multicols}{2}
\begingroup
\input{bits/varH.tex}
\endgroup
\end{multicols}
```

## HeH<sup>+</sup> 
\label{data:varHeH}

```{=latex}
\begin{multicols}{2}
\begingroup
\tiny
\input{bits/varHeH.tex}
\endgroup
\end{multicols}
```

# Código-Fonte

## Distância Variável de Ligação 
\label{code:vardist}

```python
# ---- Hydrogen molecule file ----
import psi4
import numpy as np


dmin = 0.1 # Angstrom
dmax = 1.5 # Angstrom
pt_density = 500 # Number of points in range

H = psi4.geometry('''
H
H 1 d
''')

print('# Distance (Angstrom)	Energy (Hartree)')
for r in np.linspace(dmin, dmax, pt_density):
	H.d = r
	E = psi4.energy('scf/sto-3g')
    # Print into ex1-he.out
	#print('{}	{}'.format(r, E))
```

## Ficheiros de *Input* para Espécies Ionizadas
\label{code:ioninp}

### (H2)+

```python
molecule Hion {
	1 2
	H
	H 1 0.7116232464929859
}

set basis sto-3g
set reference uhf

energy('scf', molecule=Hion)
```

### HeH(2+)

```python
molecule HeHion {
	2 2
	He
	H 1 0.92876254180602
}

set basis sto-3g
set reference uhf

energy('scf')
```

## Ficheiros Modificados para uso de UHF/cc-pVTZ
\label{code:ccpvtz}

### H2

```python
set basis cc-pVTZ 	

molecule H {
	H
	H 1 0.7116232464929859

}

energy('scf', molecule=H)
```

### HeH+

```python
molecule HeH {
	1 1
	he
	h 1 0.92876254180602 
}

set basis cc-pVTZ 	

energy('hf', molecule=HeH)
```

### (H2)+

```python
molecule Hion {
	1 2
	H
	H 1 0.7116232464929859
}

set basis cc-pvtz 
set reference uhf

energy('scf', molecule=Hion)
```

### (HeH)2+

```python
molecule HeHion {
	2 2
	He
	H 1 0.92876254180602
}

set basis cc-pvtz 
set reference uhf

energy('scf'. molecule=HeHion)
```

## Output Relevate para Cálculos com UHF/cc-pVTZ
\label{output:ccpvtz}

### H2

-------------------------------

```
    Orbital Energies [Eh]
    ---------------------

    Doubly Occupied:                                                      

       1Ag    -0.603595  

    Virtual:                                                              

       1B1u    0.171144     2Ag     0.302368     2B1u    0.654855  
       1B3u    0.695940     1B2u    0.695940     3Ag     1.104009  
       1B3g    1.124561     1B2g    1.124561     3B1u    1.538832  
       4B1u    2.343157     4Ag     2.391601     5Ag     3.180439  
       1B1g    3.180439     2B2u    3.429113     2B3u    3.429113  
       6Ag     3.684364     2B2g    3.867292     2B3g    3.867292  
       5B1u    4.116223     1Au     4.116223     3B3u    4.142879  
       3B2u    4.142879     6B1u    4.754864     3B2g    5.847560  
       3B3g    5.847560     7Ag     5.931725     7B1u    7.137166  

    Final Occupation by Irrep:
             Ag   B1g   B2g   B3g    Au   B1u   B2u   B3u 
    DOCC [     1,    0,    0,    0,    0,    0,    0,    0 ]

  Energy converged.

  @DF-RHF Final Energy:    -1.13259931280365
```

-------------------------------

### HeH+

-------------------------------

```
   Orbital Energies [Eh]
    ---------------------

    Doubly Occupied:

       1A1    -1.528449

    Virtual:

       2A1    -0.326117     3A1     0.015863     1B1     0.292834
       1B2     0.292834     4A1     0.335638     5A1     0.639348
       2B2     1.215665     2B1     1.215665     6A1     1.699397
       7A1     2.049225     1A2     2.874929     8A1     2.874929
       3B2     3.159418     3B1     3.159418     4B2     3.320952
       4B1     3.320952     9A1     3.325547    10A1     4.041023
      11A1     4.704106     2A2     5.896270    12A1     5.896270
       5B2     6.678660     5B1     6.678660     6B1     7.480816
       6B2     7.480816    13A1     7.905967    14A1     8.256953

    Final Occupation by Irrep:
             A1    A2    B1    B2
    DOCC [     1,    0,    0,    0 ]

  Energy converged.

  @DF-RHF Final Energy:    -2.92325206606155

   => Energetics <=

    Nuclear Repulsion Energy =              1.1395317635463449
    One-Electron Energy =                  -5.0686683289632297
    Two-Electron Energy =                   1.0058844993553369
    Total Energy =                         -2.9232520660615480
```
-------------------------------

### (H2)+

-------------------------------

```
    Alpha Occupied:

       1Ag    -1.304132

    Alpha Virtual:

       1B1u   -0.191396     2Ag    -0.050315     1B3u    0.193334
       1B2u    0.193334     2B1u    0.204653     3Ag     0.649944
       1B3g    0.681996     1B2g    0.681996     3B1u    1.139861
       4B1u    1.749285     4Ag     1.760757     5Ag     2.575004
       1B1g    2.575004     2B2u    2.755364     2B3u    2.755364
       6Ag     3.077701     2B2g    3.297282     2B3g    3.297282
       5B1u    3.539736     1Au     3.539736     3B3u    3.556795
       3B2u    3.556795     6B1u    4.177506     3B2g    5.222160
       3B3g    5.222160     7Ag     5.300924     7B1u    6.434774

    Beta Occupied:



    Beta Virtual:

       1Ag    -0.549582     1B1u   -0.133099     2Ag     0.006571
       1B3u    0.273697     1B2u    0.273697     2B1u    0.286977
       3Ag     0.691599     1B3g    0.705288     1B2g    0.705288
       3B1u    1.148189     4B1u    1.848159     4Ag     1.879075
       5Ag     2.626226     1B1g    2.626226     2B2u    2.859230
       2B3u    2.859230     6Ag     3.147287     2B2g    3.340278
       2B3g    3.340278     5B1u    3.568494     1Au     3.568494
       3B3u    3.590819     3B2u    3.590819     6B1u    4.210668
       3B2g    5.264352     3B3g    5.264352     7Ag     5.350657
       7B1u    6.510069

    Final Occupation by Irrep:
             Ag   B1g   B2g   B3g    Au   B1u   B2u   B3u
    DOCC [     0,    0,    0,    0,    0,    0,    0,    0 ]
    SOCC [     1,    0,    0,    0,    0,    0,    0,    0 ]

  Energy converged.

  @DF-UHF Final Energy:    -0.56051248353801
```

-------------------------------

### (HeH)2+

-------------------------------

```
    Alpha Occupied:

       1A1    -2.587460

    Alpha Virtual:

       2A1    -0.855674     3A1    -0.359847     1B2    -0.214605
       1B1    -0.214605     4A1    -0.128270     5A1     0.151636
       2B2     0.543632     2B1     0.543632     6A1     1.101810
       7A1     1.449587     1A2     2.329920     8A1     2.329920
       3B2     2.509545     3B1     2.509545     9A1     2.729528
       4B2     2.802953     4B1     2.802953    10A1     3.274551
      11A1     3.896984     2A2     5.022259    12A1     5.022259
       5B2     5.789719     5B1     5.789719     6B1     6.491647
       6B2     6.491647    13A1     6.922945    14A1     7.388450

    Beta Occupied:



    Beta Virtual:

       1A1    -1.437882     2A1    -0.763385     3A1    -0.332978
       1B1    -0.177065     1B2    -0.177065     4A1    -0.058080
       5A1     0.193199     2B1     0.632576     2B2     0.632576
       6A1     1.167996     7A1     1.485753     1A2     2.336959
       8A1     2.336959     3B2     2.535749     3B1     2.535749
       9A1     2.778168     4B2     2.808221     4B1     2.808221
      10A1     3.413885    11A1     4.004692     2A2     5.087402
      12A1     5.087402     5B2     5.871021     5B1     5.871021
       6B1     6.646467     6B2     6.646467    13A1     7.069140
      14A1     7.463575

    Final Occupation by Irrep:
             A1    A2    B1    B2
    DOCC [     0,    0,    0,    0 ]
    SOCC [     1,    0,    0,    0 ]

  Energy converged.

  @DF-UHF Final Energy:    -1.44792839770044
```

-------------------------------
