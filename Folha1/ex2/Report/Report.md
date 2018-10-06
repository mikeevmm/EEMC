---
title: Folha 1, Ex. 2 --- Estutura Electrónica e Modelação Computacional
author: Miguel Murça, 2015235874
lang: pt-PT
documentclass: article
header-includes: |
	\usepackage{supertabular, multicol, listings, physics}
---

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
\lstset{frame=single}
\begin{figure}
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
\begin{figure}
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
\begin{figure}
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
\begin{figure}
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

Pelo que, comparando a energia de ionização obtida pelo teorema de Koopmans ($I_p$) com a
energia de ionização prevista por aplicação direta de UHF/STO-3G ($I_d$), assim como 
com energias de ionização experimentalmente previstas[^CCCBDB] ($I_e$)
(tabela \ref{tab:ionenergies}): 

|            	| $H_2$     	| $HeH^+$   	| $H_2^+$ 	| $HeH^{2+}$ 	|
|------------	|-----------	|-----------	|---------	|------------	|
| E          	| -1.117553 	| -2.854418 	| -0.527029 | -1.356665     |
| $\epsilon$ 	| -0.590501    	| -1.270649 	| ---     	| ---        	|
| $I_p$      	|  0.590501   	|  1.270649  	| ---     	| ---        	|
| $I_d$      	|  0.590524  	|  1.497753  	| ---     	| ---        	|
| $I_e$         |  0.566895     | Sem Dados     | ---       | ---           |

Table: Comparação das energias HOMO ($\epsilon$), total (E), e energias
    de ionização pelo teorema de Koopmans ($I_p$) e por cálculo da energia
    da espécie ionizada ($I_d$). Todos os valores estão expressos em Hartree.
    \label{tab:ionenergies}

[^CCCBDB]:  NIST Computational Chemistry Comparison and Benchmark Database,
            NIST Standard Reference Database Number 101
            Release 19, April 2018, Editor: Russell D. Johnson III
            http://cccbdb.nist.gov/  


Observa-se que para o $H_2$ as energias de ionização $I_p$ e $I_d$ concordam em
quatro casas decimais, pelo que a energia de relaxação electrónica pode, de facto,
ser desprezada em boa aproximação. Isto é esperado, uma vez que o número de electrões
na espécie original ($H_2$) é muito reduzido, pelo que a energia de relaxação é pequena.

No entanto, observa-se que tanto $I_p$ como $I_d$ apenas concordam com $I_e$  

----------------------------------------------------------------------------------------

\twocolumn
\appendix

# *Raw Data*

## H<sub>2</sub> 
\label{data:varH}

```{=latex}
\begingroup
\input{bits/varH.tex}
\endgroup
```

## HeH<sup>+</sup> 
\label{data:varHeH}

```{=latex}
\begingroup
\tiny
\input{bits/varHeH.tex}
\endgroup
```

\onecolumn
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