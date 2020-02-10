## Volume computation and sampling

**VolEsti** is a C++ library for volume approximation and sampling of convex bodies with an *R* interface.  
#### **This is a development branch that contains the supplementary code for paper "Sampling  the feasible set of SDP and volume approximation" submitted to ISSAC 2020.**   

### - R Interface
------------

####  Install Rcpp package  
 
* Install package-dependencies: `Rcpp`, `RcppEigen`, `BH`.  

1. Then use devtools package to install `volesti` Rcpp package. In folder `/root/R-prog` Run:
```r
Rcpp::compileAttributes()  
library(devtools)  
devtools::build()  
devtools::install()  
library(volesti)  
```
2. You can use Rstudio as well to open `volesti.Rproj` and then click `build source Package` and then `Install and Restart` in `Build` at the menu bar.  


####  Run the code from `R`  
* Generate a spectrahedron using the function `generator_sdp(n,m)`. Inputs:  
  - `n` is the dimension the spectrahedron lies.  
  - `m` is the dimension of the matrices in LMI.  
  - A txt file with name `sdp_prob_n_m.txt` will be created in folder `/root/R-prog`. You cas use this file (or any other with the same format) to give it as input in the following functions.  
* Compute the volume of a spectrahedron using the function `volume()`. Inputs: 
  - `filename` is a string with the name of file in the format that the function `generator_sdp()` generates.  
  - The function calls the algorithm described in the paper to compute the volume of the input spectrahedron.  
* Sample points from a spectrahedron using the function `sample_points()`. Inputs:  
  - `file` is a string with the name of file in the format that the function `generator_sdp()` generates.  
  - `distribution` is a string to declare from which distribution to sample from: a) `uniform` or b) `boltzmann`. The default value is `uniform`.  
  - `N` is an integer to declare how many points to sample from the spectrahedron. The default value is `100`.  
  - `walk_length` is an integer to declare the walk length of the random walk. The default value is `1`.  
  - `Temperature` is a numeric input to declare the variance of the Boltzamann distribution. The default value is `1`.  
  - `random_walk` is a string to declare the random walk to use: a) `billiard` for billiard walk, b) `RDHR` for random directions Hit and Run, c) `CDHR` for coordinate directions Hit and Run or d) `HMC` for Hamiltonian Monte Carlo for reflections. The default value `billiard` for the uniform distribution and `HMC` for the Boltzmann distribution.  
The function returns a `nxN` matrix that contains the sampled points columnwise.   
* Approximate the solution of an sdp using the function `sdp_approx()`. Inputs:  
  - `filename` is a string with the name of file in the format that the function `generator_sdp()` generates.  
  - `N`is an integer to declare how many iterations to perform. The default value is `20`.  
  - `random_walk` is a string to declare the random walk to use: a) `HMC` for Hamiltonian Monte Carlo for reflections or b) `RDHR` for random directions Hit and Run. The default value is `HMC`.  
  - `walk_length` is an integer to declare the walk length of the random walk. The default value is `1`.  
The function returns a `N`-dimensional vector with the values of the objective function of each iteration.  

### - C++ Interface
------------

####  Compile C++ sources and run tests 

To compile the C++ code run in folder test:  
```
cmake .  
make  
```

