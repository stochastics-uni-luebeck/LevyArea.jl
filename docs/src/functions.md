# Function reference

```@meta
DocTestSetup = :(using LinearAlgebra; using IteratedIntegrals)
```

## Main Functions

```@docs
simiterintegrals
IteratedIntegrals.levyarea
IteratedIntegrals.AbstractIteratedIntegralAlgorithm
IteratedIntegrals.AbstractErrorNorm
optimal_algorithm
```

## Algorithmic properties

```@docs
terms_needed
IteratedIntegrals.convorder
IteratedIntegrals.errcoeff
IteratedIntegrals.norv
IteratedIntegrals.effective_cost
```

## SDE schemes

```@docs
IteratedIntegrals.em
IteratedIntegrals.sri
```

## Helper functions

```@docs
IteratedIntegrals.ito_correction!
```

## Miscellaneous

```@docs
IteratedIntegrals.createBrownianPath
IteratedIntegrals.restrict
```

## Index

```@index
```
