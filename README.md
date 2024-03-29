# MZHyper
This is a collection of MATLAB codes to reproduce all the figures and tables in our paper "On the quadrature exactness in hyperinterpolation", which is available on [arxiv 2202.13691](https://arxiv.org/abs/2202.13691).

# Interval
Please add [CVX](http://cvxr.com/cvx/), the [Chebfun](https://www.chebfun.org/) Toolbox "chebfun-master" [2], and the code files by S. Foucart for reproducing figures in [1] onto path before running the codes.

[1] COMPUTING A QUANTITY OF INTEREST FROM OBSERVATIONAL DATA
by R. DeVore, S. Foucart, G. Petrova, P. Wojtaszczyk. 
Codes are available at https://www.math.tamu.edu/~foucart/papers.html (item 31).
The zip file (reproducible+supplement) also contains CVX.

# Sphere
Please add the sphere_approx_toolbox_v3.0 (can be found in this repository) onto path before running the codes.

Please go to /sphere_approx_toolbox_v3.0/utilities/ and change the bold part of the path '**/Users/haoningwu/Documents/MATLAB/hyperinterpolation**/sphere_approx_toolbox_v3.0/data/xx' in **loadStd.m** to your own path storing the sphere_approx_toolbox_v3.0. Otherwise, MATLAB would report error:
  >The file '/Users/haoningwu/Documents/MATLAB/hyperinterpolation/sphere_approx_toolbox_v3.0/data/xx/xxxx' could not be opened because: No such file or
directory
