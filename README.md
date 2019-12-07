## Supplementary code for paper "Computing the volume of projections of polytopes" submitted to SoCG 2020

#### 0. Compile code
- Clone repository and switch to branch socg_2020:  
```
git clone https://github.com/GeomScale/volume_approximation.git  
git checkout socg_2020
```
- Save `liblpsolve55.so` in folder `/usr/lib/lp_solve/`. You will find it in `/external` folder.  
- In folder `/test` compile C++ code by running:  
```
cmake .  
make  
```

#### 1. Zonotopes (Table 1):  
- You can generate a random zonotope in dimension `dim` with `n` generators. Use flag `-dist` to select generator (the default is the uniform zonotope): `-dist` 1 for uniform, `-dist` 2 for Gaussian, `-dist` 3 for exponential. Run:  
```
./generate -zonotope -d dim -k n -dist 1
```

- Estimate the volume using balls in MMC:  
```
./vol -f3 zonotope_dim_k.ext
```

- Estimate the volume using h-polytopes in MMC:  
```
./vol -f3 zonotope_dim_k.ext -hpoly
```

- For example the following commands:  
```
./generate -zonotope -d 10 -m 15
./vol -f3 zonotope_10_15.ext -hpoly
```
Will generate a uniform 10-dimensional zonotope with 15 generators and estimate the volume by using h-polytopes in MMC.  

- To compute the exact volume run:  
```
./vol -f3 zonotope_10_15.ext -exact_zono
```

#### 2. V-polytopes (Table 3):  

- To generate a cross polytope in `dim` dimension and estimate the volume:  
```
./generate -cross -v -d dim
./vol -f2 cross_dim.ext
```

- To generate a unit simplex in `dim` dimension:  
```
./generate -simplex -v -d dim
./vol -f2 simplex_dim.ext
```

- To generate a unit cube in `dim` dimension and estimate the volume:  
```
./generate -cube -v -d dim
./vol -f2 cube_dim.ext
```

- Generate a random V-polytope in `dim` dimension with `n` vertices. Use flag `-r` for rounding and the flag `-body` to select generator (the default is `rvs`. `-body` 1 is for `rvs` and `-body` 2 is for `rvc`. Estimate the volume:  
```
./generate -v -rand -body 1 -d dim -k n
./vol -f2 rvs_dim_k.ext -r
```

- Estimate the volume using h-polytopes in MMC:  
```
./vol -f2 rvs_dim_k.ext -hpoly
```

Note: If you wish to give a specific V- or Z-polytope as input use an `.ext` file. Keep the same format as in the generated files.

#### 1. H-polytopes:  
- To generate a unit cube in `dim` dimension and estimate the volume:  
```
./generate -cube -h -d dim
./vol -f1 cube_dim.ine
```
For example:  
```
./generate -cube -h -d 20
./vol -f1 cube_20.ine
```

- To generate a random H-polytope in dimension `dim` with `n` facets. Again use flag `-r` for rounding. Estimate the volume:  
```
./generate -h -rand -d dim -k n
./vol -f1 rhs_dim_k.ine
```


#### 4. Flags

- Use flag `-rdhr` to use Random Directions HnR:  
```
./vol -f3 zonotope_dim_k.ext -hpoly -rdhr
```

- Use flag `-WalkL` to set the step of the random walk (the default value is 1). For example:  
```
./vol -f3 zonotope_dim_k.ext -WalkL 5
```
Will set the step equals to 5.  


- Use flag `-e` to set the error (the default value is `0.1`). For example:  
```
./vol -f3 zonotope_dim_k.ext -e 0.2
```

- Use flag `-WinLen` to set the length of the sliding window (the default value is 125). For example:  
```
./vol -f2 random_v_poly_dim_k.ext -r -WinLen 500
```
Will set the window's length `l=500`.  

- Use flags `-l` and `-u` to set the test values for L-test (r) and U-test (r + \delta) respectively. For example:  
```
./vol -f3 zonotope_dim_k.ext -l 0.01 -u 0.05
```
Will define ratios between `0.01` and `0.015` with high probability.  

- Use flag `-nuN` to set the number of points that are generated in each step of the annealing schedule, from the convex body P_i of the previous step. For example:  
```
./vol -f3 zonotope_dim_k.ext -nuN 1600 10
```
Wil sample 1600 points in total and split them to 10 sub-lists. So the degrees of freedom in each t-test will be 9 = 10-1.  

#### 5. Test PCA over-aproximations of a zonotope

- To compute the ratio for the PCA approximation of a zonotope that is described in a `.ext` file, use flag `-pca` and run:  
```
./vol -f3 zonotope_dim_k.ext -hpoly -pca
```

