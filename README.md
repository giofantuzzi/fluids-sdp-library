# Fluid mechanics SDP library
This repository is a library of SDPs obtained from a problem in fluid mechanics. Background on the problem can be found at https://quinopt.readthedocs.io/04_examples/planeCouetteBF.html.

**NOTE:** The SDPs in this library are **not** meant to give accurate solutions of the fluid flow problem described at the above link. Rather, they are a testbed for low-rank SDP solvers and warm-starting strategies.

## Problem & file description
The SDPs in this library are indexed by a parameter Re (the Reynolds number), which increases in logarithmically spaced steps from 100 to 1000. The optimal solution varies continuously as Re is increased.

At each value of Re, the optimal positive semidefinite matrices in the primal form problem
```
    min  c'*x
    s.t. A*x = b, x \in K
```
have rank at most one.

## Challenges
1. Speed up the solution of each SDP by exploiting the low-rank structure of the primal-standard-form optimizers
2. Effectively solve the whole sequence of SDPs using warm-starts 

## Contents of each file
Each .mat file in this reposity contains the following variables.

1. `Re`: The Reynolds number (parameter) for which the SDP was generated.
2. `SDP`: A structure containing the SDP in SeDuMi format. The fields `SDP.A`, `SDP.b`, `SDP.c`, `SDP.K` specify the standard forms of a conic program as stated above. The only cone in these programs is a direct product of semidefinite cones whose size is listed in `SDP.K.s`.