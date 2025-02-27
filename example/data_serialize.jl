using TSSOS, DynamicPolynomials
using Serialization
using JuMP

@polyvar x[1:3]

f1 = x[1]^4 + (x[1] * x[2] - 1)^2
f2 = x[2]^2 * x[3]^2 + (x[3]^2 - 1)^2

f = f1 + f2

opt, sol, data = cs_tssos_first([f], x, 2, CS=false, TS=false, solver= "COSMO")

cs_tssos_first_outer_dense = deserialize("./example/data/nblockmix/cs_tssos_first_outer.jls")