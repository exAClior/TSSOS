using TSSOS, DynamicPolynomials
using Serialization
using JuMP

@polyvar x[1:3]

f1 = x[1]^4 + (x[1] * x[2] - 1)^2
f2 = x[2]^2 * x[3]^2 + (x[3]^2 - 1)^2

f = f1 + f2
arch = 3 - sum(x[1]^2 + x[2]^2)
arch2 = 3 - sum(x[2]^2 + x[3]^2)

opt, sol, data = cs_tssos_first([f], x, 2, CS=false, TS=false, solver= "COSMO")
# opt, sol, data = cs_tssos_first([f;arch;arch2], x, 2, CS=false, TS=false, solver= "COSMO")

# need to use COSMO otherwise this will stop due to slow progress
# opt, sol, data = cs_tssos_first([f;arch;arch2], x, 2, TS=false, solver= "COSMO")
opt, sol, data = cs_tssos_first([f], x, 2, TS=false, solver= "COSMO")

# load data needed to construct sdp
dense_data = deserialize("./example/data/dense_solvesdp_input_data.jls")
cs_data = deserialize("./example/data/cs_solvesdp_input_data.jls")

# remove all normality = true branch
@assert !dense_data["normality"]
@assert !cs_data["normality"]

# remove all MomentOne = true || solution = true branch
@assert !dense_data["MomentOne"]
@assert !cs_data["MomentOne"]
@assert !dense_data["solution"]
@assert !cs_data["solution"]

# remove all TS = true branch
@assert !dense_data["TS"]
@assert !cs_data["TS"]

@assert  !(dense_data["numeq"] > 0)
@assert  !(cs_data["numeq"] > 0)

@assert !dense_data["tune"]
@assert !cs_data["tune"]

@assert !dense_data["dualize"]
@assert !cs_data["dualize"]

@assert dense_data["solve"]
@assert cs_data["solve"]


dense_data["blocksize"]
cs_data["blocksize"]

dense_data["cql"]
cs_data["cql"]

dense_data["cl"]
cs_data["cl"]

dense_data["cliquesize"]
cs_data["cliquesize"]

opt, ksupp, moment_dense, GramMat_dense, multiplier_equality, SDP_status, model_dense = TSSOS.solvesdp_debug(dense_data["n"], dense_data["m"], dense_data["supp"], dense_data["coe"], dense_data["basis"], dense_data["hbasis"], dense_data["cliques"], dense_data["cql"], dense_data["cliquesize"], dense_data["I"], dense_data["J"], dense_data["ncc"], dense_data["blocks"], dense_data["eblocks"], dense_data["cl"], dense_data["blocksize"], numeq=dense_data["numeq"], nb=dense_data["nb"], QUIET=false, signsymmetry=dense_data["signsymmetry"], TS=dense_data["TS"], solver="COSMO", tune=dense_data["tune"], dualize=dense_data["dualize"], solve=dense_data["solve"], solution=true, MomentOne=dense_data["MomentOne"], Gram=true, Mommat=true, cosmo_setting=dense_data["cosmo_setting"], mosek_setting=dense_data["mosek_setting"], normality=dense_data["normality"], NormalSparse=dense_data["NormalSparse"])

@assert is_solved_and_feasible(model_dense)

opt, ksupp, moment_cs, GramMat_cs, multiplier_equality, SDP_status, model_cs = TSSOS.solvesdp_debug(cs_data["n"], cs_data["m"], cs_data["supp"], cs_data["coe"], cs_data["basis"], cs_data["hbasis"], cs_data["cliques"], cs_data["cql"], cs_data["cliquesize"], cs_data["I"], cs_data["J"], cs_data["ncc"], cs_data["blocks"], cs_data["eblocks"], cs_data["cl"], cs_data["blocksize"], numeq=cs_data["numeq"], nb=cs_data["nb"], QUIET=false, signsymmetry=cs_data["signsymmetry"], TS=cs_data["TS"], solver="COSMO", tune=cs_data["tune"], dualize=cs_data["dualize"], solve=cs_data["solve"], solution=true, MomentOne=cs_data["MomentOne"], Gram=true, Mommat=true, cosmo_setting=cs_data["cosmo_setting"], mosek_setting=cs_data["mosek_setting"], normality=cs_data["normality"], NormalSparse=cs_data["NormalSparse"])

@assert is_solved_and_feasible(model_cs)

objective_value(model_cs)



cs_data["cliques"]

cs_clique1_vars = x[cs_data["cliques"][1]]
cs_clique2_vars = x[cs_data["cliques"][2]]

cs_clique1_basis = vcat([monomials(cs_clique1_vars, max_degree) for max_degree in 0:2]...)

cs_clique1_monomial2moment_matrix_idx = Dict([mono_r * mono_c => (i, j) for (i, mono_r) in enumerate(cs_clique1_basis), (j, mono_c) in enumerate(cs_clique1_basis)])

cs_clique2_basis = vcat([monomials(cs_clique2_vars, max_degree) for max_degree in 0:2]...)

cs_clique2_monomial2moment_matrix_idx = Dict([mono_r * mono_c => (i, j) for (i, mono_r) in enumerate(cs_clique2_basis), (j, mono_c) in enumerate(cs_clique2_basis)])

reconstructed_moment_cs_total = zeros(size(moment_dense[1]))

dense_basis = vcat([monomials(x, max_degree) for max_degree in 0:2]...)
for (i,mono_r) in enumerate(dense_basis), (j,mono_c) in enumerate(dense_basis)
	if haskey(cs_clique1_monomial2moment_matrix_idx, mono_r * mono_c)
		reconstructed_moment_cs_total[i,j] += moment_cs[1][cs_clique1_monomial2moment_matrix_idx[mono_r * mono_c]...]
	end
	if haskey(cs_clique2_monomial2moment_matrix_idx, mono_r * mono_c)
		reconstructed_moment_cs_total[i,j] += moment_cs[2][cs_clique2_monomial2moment_matrix_idx[mono_r * mono_c]...]
	end
end

using LinearAlgebra

reconstructed_moment_cs_total

eigvals(reconstructed_moment_cs_total)

moment_dense[1]





display(round.(moment_cs[1], digits=10))
display(round.(moment_cs[2], digits=10))

display(round.(moment_dense[1], digits=10))

# corresponds to variable pos


dense_data["cliques"]
cs_data["cliques"]

display(round.(GramMat_dense[1][1][1], digits=8))
display(round.(GramMat_cs[1][1][1], digits=8))
display(round.(GramMat_cs[2][1][1], digits=8))


# opt, sol, data = tssos_first([f], x, 2, TS="block", solver="COSMO")
# opt, sol, data = tssos_higher!(data, TS="block")