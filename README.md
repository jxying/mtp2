# mtp2 
Matlab implementation of fast projected Newton-like (FPN) method [1] for learning MTP2 Gaussian graphical models. The problem can be formulated as

$$
\underset{\boldsymbol{\Theta}}{\mathsf{minimize}}  -\log\det\left(\boldsymbol{\Theta}\right)+\left\langle \boldsymbol{\Theta},\mathbf{S}\right\rangle +\sum_{i\neq j}\Lambda_{ij}\left|\Theta_{ij}\right|, 
$$

subject to  

$$ 
	\boldsymbol{\Theta}\succ\mathbf{0}, \text{ and } \Theta_{ij}\leq0,\forall i\neq j
$$ 

## To run the code:
(1) Download the source files.
(2) Run 'demo.m'

You may consider using the FPN solver with [bridge-block decomposition](https://github.com/jxying/mtp2-bbd) for learning large-scale MTP2 Gaussian graphical models. The bridge-block decomposition is designed to reduce the computational and memory costs of existing algorithms like FPN, particularly in the context of large-scale data.

## References

[1] Jian-Feng Cai, José Vinícius de Miranda Cardoso, Daniel P. Palomar, and Jiaxi Ying, "Fast Projected Newton-like Method for Precision Matrix Estimation under Total Positivity", Advances in Neural Information Processing Systems (NeurIPS), New Orleans, LA, USA, Dec. 2023.
