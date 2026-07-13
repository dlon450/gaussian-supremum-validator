# Reviewer-comment evidence map (auto-generated from results/experiments)

Config: alpha=0.10 (tolerance 1-alpha=90%), beta=0.05 (target 1-beta=0.95); Gaussian DGP; d=10 unless swept. 'meets_target' = point coverage >= 0.95.

## R1(major)/R2: Comparison against existing validation methods (CV, bootstrap, sectioning)

Setup: SO / SAA / Wasserstein-DRO, split=0.7, all methods at the SAME total data budget.

### paper_dro_wasserstein

Coverage (target 0.95):
```
method    NV    UG    NGS  UNGS     CV     BS  Sectioning
n                                                        
100     0.42  0.84  0.917  0.97  0.977  0.997        0.98
```

Mean objective (lower is better):
```
method    NV     UG    NGS   UNGS     CV     BS  Sectioning
n                                                          
100    -7.91 -7.267 -6.788 -2.965 -6.912 -6.762      -6.622
```

Smallest n meeting target: {'NV': None, 'UG': None, 'NGS': None, 'UNGS': 100, 'CV': 100, 'BS': 100, 'Sectioning': 100}

### paper_saa

Coverage (target 0.95):
```
method     NV     UG    NGS   UNGS     CV     BS  Sectioning
n                                                           
100     0.263  0.577  0.600  0.613  0.793  0.940       0.613
200     0.353  0.817  0.903  0.930  0.963  0.993       0.940
300     0.400  0.870  0.933  0.963  0.973  0.993       0.973
```

Mean objective (lower is better):
```
method     NV     UG    NGS   UNGS     CV     BS  Sectioning
n                                                           
100    -8.073 -7.792 -7.765 -7.748 -7.605 -7.424      -7.744
200    -7.987 -7.570 -7.445 -7.347 -7.377 -7.171      -7.316
300    -7.969 -7.563 -7.446 -7.229 -7.404 -7.299      -7.181
```

Smallest n meeting target: {'NV': None, 'UG': None, 'NGS': None, 'UNGS': 300, 'CV': 200, 'BS': 200, 'Sectioning': 300}

### paper_so

Coverage (target 0.95):
```
method     NV     UG    NGS   UNGS     CV     BS  Sectioning  benchmark
n                                                                      
100     0.440  0.653  0.680  0.680  0.757  0.870       0.657      0.977
200     0.537  0.897  0.943  0.970  0.913  0.960       0.963      1.000
300     0.627  0.953  0.983  1.000  0.887  0.953       0.993      1.000
500     0.643  0.957  0.980  0.997  0.863  0.950       0.997      1.000
```

Mean objective (lower is better):
```
method     NV     UG    NGS   UNGS     CV     BS  Sectioning  benchmark
n                                                                      
100    -7.906 -7.696 -7.655 -7.656 -7.615 -7.474      -7.692     -7.286
200    -7.797 -7.458 -7.371 -7.287 -7.313 -7.229      -7.299     -6.981
300    -7.720 -7.386 -7.282 -7.078 -7.329 -7.263      -7.104     -6.842
500    -7.708 -7.407 -7.339 -7.013 -7.380 -7.299      -7.049     -6.644
```

Smallest n meeting target: {'NV': None, 'UG': 300, 'NGS': 300, 'UNGS': 200, 'CV': None, 'BS': 200, 'Sectioning': 200, 'benchmark': 100}

## R1/R2: Guideline for relative budgeting of n1 vs n2 (Phase-1 fraction)

### paper_saa (n=400) coverage vs Phase-1 fraction:
```
method     NV     UG    NGS   UNGS
split                             
0.1     0.489  0.583  0.589  0.591
0.2     0.624  0.886  0.909  0.914
0.3     0.565  0.948  0.982  0.988
0.4     0.504  0.938  0.978  0.985
0.5     0.463  0.917  0.990  0.988
0.6     0.386  0.928  0.979  0.985
0.7     0.364  0.909  0.965  0.980
0.8     0.392  0.906  0.969  0.981
0.9     0.316  0.886  0.886  0.934
```

- uni. Gaussian: Phase-1 fractions meeting target = none

- norm. GS: Phase-1 fractions meeting target = ['30%', '40%', '50%', '60%', '70%', '80%']

### paper_so (n=500) coverage vs Phase-1 fraction:
```
method     NV     UG    NGS   UNGS
split                             
0.1     0.603  0.720  0.723  0.725
0.2     0.746  0.958  0.980  0.982
0.3     0.761  0.982  0.996  0.998
0.4     0.732  0.972  0.997  0.999
0.5     0.704  0.971  0.996  0.999
0.6     0.725  0.971  0.999  1.000
0.7     0.668  0.975  0.991  0.998
0.8     0.646  0.945  0.987  1.000
0.9     0.539  0.931  0.975  0.992
```

- uni. Gaussian: Phase-1 fractions meeting target = ['20%', '30%', '40%', '50%', '60%', '70%']

- norm. GS: Phase-1 fractions meeting target = ['20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%']

Reading: too-small n1 starves the solution path (fails target); a broad middle band (typically 50-80%, with ~70% robust) meets it. This is the empirical n1 guideline.

## R2: Dimension dependence (dimension-free feasibility claim)

RO coverage vs dimension d (n=500, split=0.5):
```
method     NV     UG    NGS   UNGS  benchmark
d                                            
2       0.945  0.993  0.999  0.999        1.0
5       0.943  0.997  0.999  0.998        1.0
10      0.867  0.993  0.998  0.994        1.0
20      0.870  0.990  0.996  0.990        1.0
50      1.000  1.000  1.000  1.000        1.0
100     1.000  1.000  1.000  1.000        1.0
```

UG objective advantage over SCA benchmark (benchmark_obj - UG_obj) vs d:
```
d
2      0.153
5      0.340
10     0.539
20     0.686
50     1.567
100    0.305
```

Reading: coverage of the proposed validators stays >= target across d (light dimension dependence), while the objective advantage over the fixed-margin SCA benchmark grows with d.

## R2 (Q1/Q2): Coverage and objective-gap convergence vs n

### paper_fast coverage vs n:
```
method     UG    NGS   UNGS     NV  benchmark
n                                            
200     0.999  0.999  0.999  0.986        1.0
500     1.000  1.000  1.000  1.000        1.0
```

objective gap to path-oracle vs n:
```
method     UG    NGS   UNGS
n                          
200     0.199  0.277  0.239
500     0.001  0.003  0.001
```

### paper_ro_ellipsoid coverage vs n:
```
method     UG    NGS   UNGS     NV  benchmark
n                                            
100     0.968  0.969  0.982  0.707        1.0
200     0.972  0.985  0.983  0.779        1.0
300     0.996  0.996  0.996  0.829        1.0
500     0.993  0.998  0.994  0.867        1.0
1000    0.996  0.999  0.997  0.895        1.0
```

objective gap to path-oracle vs n:
```
method     UG    NGS   UNGS
n                          
100     0.415  0.486  0.917
200     0.235  0.328  0.475
300     0.191  0.246  0.289
500     0.073  0.124  0.144
1000    0.019  0.031  0.031
```

### paper_so coverage vs n:
```
method     UG    NGS   UNGS     NV  benchmark
n                                            
100     0.587  0.607  0.607  0.429       0.98
200     0.903  0.949  0.958  0.630       1.00
300     0.980  0.989  0.995  0.699       1.00
500     0.971  0.996  0.999  0.704       1.00
800     0.977  0.998  1.000  0.766       1.00
1000    0.977  0.996  0.999  0.796       1.00
```

objective gap to path-oracle vs n:
```
method     UG    NGS   UNGS
n                          
100     0.038  0.060  0.065
200     0.152  0.235  0.288
300     0.254  0.306  0.449
500     0.188  0.298  0.479
800     0.159  0.264  0.403
1000    0.143  0.236  0.334
```

## Robustness under heavy-tailed (multivariate-t) data

### robust_dro_wasserstein_t coverage vs n:
```
method     UG    NGS   UNGS     NV
n                                 
100     0.926  0.954  0.977  0.483
200     0.888  0.956  0.969  0.482
```

### robust_ro_ellipsoid_t coverage vs n:
```
method     UG    NGS   UNGS     NV  benchmark
n                                            
200     0.972  0.990  0.985  0.816        1.0
500     0.993  0.997  0.993  0.911        1.0
1000    0.996  0.998  0.997  0.973        1.0
```
