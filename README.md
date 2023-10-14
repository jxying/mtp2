# mtp2-bbd-Pypkg 
Matlab implementation of fast projected Newton-like (FPN) method [1] for learning large-scale MTP2 Gaussian graphical models and its accleration by introducing bridge-block decomposition [2]. The problem can be formulated as

$$
\mathsf{minimize}  -\log\det\left(\boldsymbol{\Theta}\right)+\left\langle \boldsymbol{\Theta},\mathbf{S}\right\rangle +\sum_{i\neq j}\Lambda_{ij}\left|\Theta_{ij}\right|, 
$$

subject to  

$$ 
	\boldsymbol{\Theta}\succ\mathbf{0}, \text{ and } \Theta_{ij}\leq0,\forall i\neq j
$$ 

To run the code:
(1) Download the source files.
(2) Run 'demo.m'
