# Comprehensive experiments — analysis

Config: alpha=0.10 (tolerance 90%), beta=0.05 (target coverage 0.95); Gaussian DGP.

## Best Phase-1 fraction n1/n vs N and D (budget matrix)

For each (formulation, N, D) and validator, the FEASIBLE BAND is the set of n1 fractions with coverage >= 0.95; the RECOMMENDED n1 is the feasible fraction with the lowest (best) objective (ties -> smallest n1).


### Recommended n1 (best feasible), UG

```
n               100   200   400   500   1000
formulation d                               
RO          2    NaN   0.3   NaN   0.2   0.2
            10   NaN   0.1   NaN   0.2   0.1
            50   NaN   0.1   NaN   0.1   0.1
SAA         5    NaN   NaN   0.3   NaN   NaN
SO          2    NaN   0.3   NaN   0.1   0.1
            10   NaN   NaN   NaN   0.2   0.1
            50   NaN   NaN   NaN   0.3   0.2
Wasserstein 10   0.3   0.3   NaN   NaN   NaN
            5    0.2   0.2   NaN   NaN   NaN
```


### Recommended n1 (best feasible), NGS

```
n               100   200   400   500   1000
formulation d                               
RO          2    NaN   0.2   NaN   0.3   0.1
            10   NaN   0.2   NaN   0.2   0.2
            50   NaN   0.1   NaN   0.1   0.2
SAA         5    NaN   0.5   0.4   NaN   NaN
            10   NaN   0.6   0.3   NaN   NaN
            20   NaN   NaN   0.3   NaN   NaN
SO          2    NaN   0.3   NaN   0.1   0.1
            10   NaN   0.5   NaN   0.2   0.1
            50   NaN   NaN   NaN   0.3   0.2
Wasserstein 5    0.4   0.4   NaN   NaN   NaN
            10   0.3   0.4   NaN   NaN   NaN
```

Reading: entries are the recommended Phase-1 fraction; `NaN` = no tested n1 met target (data too small). Trends to note: how the recommended fraction moves as N grows (more data -> smaller n1 can suffice) and as D grows (higher dim -> need more Phase-1 data to estimate the path).


## Feasibility & optimality vs N and D (ndgrid matrix, split=0.5)


### UG coverage over (N x D)

Wasserstein:
```
n     100    200    400
d                      
2   0.945  0.928  0.958
5   0.899  0.914  0.923
10  0.918  0.928  0.928
20  0.932  0.926  0.931
```

FAST:
```
n    200   500   1000
d                    
2   1.000   1.0   1.0
5   1.000   1.0   1.0
10  0.999   1.0   1.0
20  0.995   1.0   1.0
50  0.992   1.0   1.0
```

RO:
```
n     100    200    500    1000   2000
d                                     
2    0.972  0.977  0.993  1.000  1.000
5    0.961  0.986  0.997  1.000  1.000
10   0.968  0.972  0.993  0.996  1.000
20   0.958  0.984  0.990  0.994  0.999
50   0.973  0.994  1.000  1.000  1.000
100  0.975  1.000  1.000  1.000  1.000
```

SAA:
```
n     100    200    400
d                      
2   0.894  0.932  0.934
5   0.730  0.913  0.933
10  0.600  0.877  0.917
20  0.481  0.829  0.909
```

SO:
```
n    100    200    500    1000   2000
d                                    
2   0.908  0.968  0.990  0.992  0.999
5   0.769  0.949  0.978  0.987  0.999
10  0.587  0.905  0.975  0.985  0.990
20  0.482  0.842  0.975  0.975  0.989
50  0.488  0.823  0.967  0.982  0.991
```


### NGS coverage over (N x D)

Wasserstein:
```
n     100    200    400
d                      
2   0.951  0.962  0.991
5   0.930  0.971  0.978
10  0.949  0.973  0.986
20  0.956  0.967  0.993
```

FAST:
```
n    200   500   1000
d                    
2   1.000   1.0   1.0
5   1.000   1.0   1.0
10  0.999   1.0   1.0
20  0.999   1.0   1.0
50  0.994   1.0   1.0
```

RO:
```
n     100    200    500    1000   2000
d                                     
2    0.973  0.992  0.999  1.000  1.000
5    0.963  0.991  0.999  1.000  1.000
10   0.969  0.985  0.998  0.999  1.000
20   0.958  0.994  0.996  0.998  0.999
50   0.973  0.995  1.000  1.000  1.000
100  0.975  1.000  1.000  1.000  1.000
```

SAA:
```
n     100    200    400
d                      
2   0.899  0.974  0.981
5   0.736  0.962  0.987
10  0.610  0.927  0.988
20  0.491  0.878  0.978
```

SO:
```
n    100    200    500    1000   2000
d                                    
2   0.917  0.989  0.997  0.999  1.000
5   0.786  0.980  0.988  1.000  1.000
10  0.607  0.948  0.996  0.997  0.998
20  0.496  0.885  0.994  0.997  0.998
50  0.496  0.859  0.988  0.996  0.998
```


## CV / bootstrap at K in {3,5,10} vs proposed (folds matrix, split=0.7)


### Wasserstein coverage vs N (target 0.95)
```
method    UG    NGS  UNGS    CV3    CV5  CV10    BS3   BS5   BS10   Sec3  Sec5  Sec10
n                                                                                    
100     0.84  0.917  0.97  0.943  0.937  0.98  0.933  0.98  0.997  0.647  0.98   0.96
```


### SAA coverage vs N (target 0.95)
```
method     UG    NGS   UNGS   CV3    CV5   CV10    BS3    BS5   BS10   Sec3   Sec5  Sec10
n                                                                                        
100     0.577  0.600  0.613  0.69  0.723  0.763  0.733  0.803  0.933  0.370  0.567  0.637
200     0.817  0.903  0.930  0.92  0.890  0.967  0.870  0.960  0.997  0.633  0.767  0.940
```


### SO coverage vs N (target 0.95)
```
method     UG    NGS   UNGS    CV3   CV5   CV10    BS3    BS5   BS10   Sec3   Sec5  Sec10  benchmark
n                                                                                                   
100     0.660  0.680  0.680  0.690  0.68  0.700  0.650  0.757  0.847  0.480  0.657  0.683      0.977
200     0.900  0.943  0.973  0.870  0.85  0.893  0.817  0.897  0.963  0.713  0.843  0.943      1.000
300     0.957  0.983  1.000  0.857  0.84  0.900  0.753  0.863  0.960  0.817  0.977  0.987      1.000
500     0.957  0.983  0.997  0.863  0.84  0.887  0.753  0.857  0.917  0.883  0.953  1.000      1.000
```
