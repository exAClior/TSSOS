using TSSOS, DynamicPolynomials
using Serialization
using JuMP

@polyvar x[1:3]

f1 = x[1]^4 + (x[1] * x[2] - 1)^2
f2 = x[2]^2 * x[3]^2 + (x[3]^2 - 1)^2

f = f1 + f2

opt, sol, data = cs_tssos_first([f], x, 2, CS=false, TS=false)

opt, sol, data = cs_tssos_first([f], x, 2, TS=false)

dense_data = deserialize("./example/data/dense_solvesdp_input_data.jls")

opt, ksupp, moment, GramMat, multiplier_equality, SDP_status, model1 = TSSOS.solvesdp_debug(dense_data["n"], dense_data["m"], dense_data["supp"], dense_data["coe"], dense_data["basis"], dense_data["hbasis"], dense_data["cliques"], dense_data["cql"], dense_data["cliquesize"], dense_data["I"], dense_data["J"], dense_data["ncc"], dense_data["blocks"], dense_data["eblocks"], dense_data["cl"], dense_data["blocksize"], numeq=dense_data["numeq"], nb=dense_data["nb"], QUIET=false, signsymmetry=dense_data["signsymmetry"], TS=dense_data["TS"], solver=dense_data["solver"], tune=dense_data["tune"], dualize=dense_data["dualize"], solve=dense_data["solve"], solution=dense_data["solution"], MomentOne=dense_data["MomentOne"], Gram=dense_data["Gram"], Mommat=dense_data["Mommat"], cosmo_setting=dense_data["cosmo_setting"], mosek_setting=dense_data["mosek_setting"], normality=dense_data["normality"], NormalSparse=dense_data["NormalSparse"])

model1

is_solved_and_feasible(model1)


cs_data = deserialize("./example/data/cs_solvesdp_input_data.jls")

opt, ksupp, moment, GramMat, multiplier_equality, SDP_status, model2 = TSSOS.solvesdp_debug(cs_data["n"], cs_data["m"], cs_data["supp"], cs_data["coe"], cs_data["basis"], cs_data["hbasis"], cs_data["cliques"], cs_data["cql"], cs_data["cliquesize"], cs_data["I"], cs_data["J"], cs_data["ncc"], cs_data["blocks"], cs_data["eblocks"], cs_data["cl"], cs_data["blocksize"], numeq=cs_data["numeq"], nb=cs_data["nb"], QUIET=false, signsymmetry=cs_data["signsymmetry"], TS=cs_data["TS"], solver=cs_data["solver"], tune=cs_data["tune"], dualize=cs_data["dualize"], solve=cs_data["solve"], solution=cs_data["solution"], MomentOne=cs_data["MomentOne"], Gram=cs_data["Gram"], Mommat=cs_data["Mommat"], cosmo_setting=cs_data["cosmo_setting"], mosek_setting=cs_data["mosek_setting"], normality=cs_data["normality"], NormalSparse=cs_data["NormalSparse"])

objective_value(model2)

num_variables(model2)
is_solved_and_feasible(model2)




# opt, sol, data = tssos_first([f], x, 2, TS="block")
# opt, sol, data = tssos_higher!(data, TS="block")