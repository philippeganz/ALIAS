# ASTROQUT

Astrophysicists are interested in recovering the 3D gas emissivity of a galaxy cluster from a 2D image taken by a telescope. A blurring phenomenon and presence of point sources make this inverse problem even harder to solve. The current state-of-the-art technique is two step: First identify the location of potential point sources, then mask these locations and deproject the data.  
 
We instead model the data as a Poisson generalized linear model (involving blurring, Abel and wavelets operators) regularized by two lasso penalties to induce sparse wavelet representation and sparse point sources. The amount of sparsity is controlled by two quantile universal thresholds. As a result, our method outperforms the existing one.
