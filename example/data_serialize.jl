using TSSOS, DynamicPolynomials
using Serialization
using JuMP

@polyvar x[1:3]

f1 = x[1]^4 + (x[1] * x[2] - 1)^2
f2 = x[2]^2 * x[3]^2 + (x[3]^2 - 1)^2

f = f1 + f2

opt, sol, data = cs_tssos_first([f], x, 2, CS=false, TS=false, solver= "COSMO")

cs_tssos_first_outer_dense = deserialize("./example/data/nblockmix/cs_tssos_first_outer.jls")


cs_tssos_first_outer_dense["supp_before_cs_tssos_first"]
cs_tssos_first_outer_dense["coe_before_cs_tssos_first"]
cs_tssos_first_outer_dense["n_before_cs_tssos_first"]
cs_tssos_first_outer_dense["d_before_cs_tssos_first"]
cs_tssos_first_outer_dense["numeq_before_cs_tssos_first"]
cs_tssos_first_outer_dense["nb_before_cs_tssos_first"]
cs_tssos_first_outer_dense["CS_before_cs_tssos_first"]
cs_tssos_first_outer_dense["cliques_before_cs_tssos_first"]
cs_tssos_first_outer_dense["basis_before_cs_tssos_first"]
cs_tssos_first_outer_dense["hbasis_before_cs_tssos_first"]
cs_tssos_first_outer_dense["minimize_before_cs_tssos_first"]
cs_tssos_first_outer_dense["TS_before_cs_tssos_first"]
cs_tssos_first_outer_dense["merge_before_cs_tssos_first"]
cs_tssos_first_outer_dense["md_before_cs_tssos_first"]
cs_tssos_first_outer_dense["QUIET_before_cs_tssos_first"]
cs_tssos_first_outer_dense["solver_before_cs_tssos_first"]
cs_tssos_first_outer_dense["tune_before_cs_tssos_first"]
cs_tssos_first_outer_dense["dualize_before_cs_tssos_first"]
cs_tssos_first_outer_dense["solve_before_cs_tssos_first"]
cs_tssos_first_outer_dense["solution_before_cs_tssos_first"]
cs_tssos_first_outer_dense["Gram_before_cs_tssos_first"]
cs_tssos_first_outer_dense["MomentOne_before_cs_tssos_first"]
cs_tssos_first_outer_dense["Mommat_before_cs_tssos_first"]
cs_tssos_first_outer_dense["tol_before_cs_tssos_first"]
cs_tssos_first_outer_dense["cosmo_setting_before_cs_tssos_first"]
cs_tssos_first_outer_dense["mosek_setting_before_cs_tssos_first"]
cs_tssos_first_outer_dense["normality_before_cs_tssos_first"]
cs_tssos_first_outer_dense["NormalSparse_before_cs_tssos_first"]


opt, sol, data = cs_tssos_first(
	cs_tssos_first_outer_dense["supp_before_cs_tssos_first"],
    cs_tssos_first_outer_dense["coe_before_cs_tssos_first"],
    cs_tssos_first_outer_dense["n_before_cs_tssos_first"],
    cs_tssos_first_outer_dense["d_before_cs_tssos_first"],
    numeq=cs_tssos_first_outer_dense["numeq_before_cs_tssos_first"],
    nb=cs_tssos_first_outer_dense["nb_before_cs_tssos_first"],
    CS=cs_tssos_first_outer_dense["CS_before_cs_tssos_first"],
    cliques=cs_tssos_first_outer_dense["cliques_before_cs_tssos_first"],
    basis=cs_tssos_first_outer_dense["basis_before_cs_tssos_first"],
    hbasis=cs_tssos_first_outer_dense["hbasis_before_cs_tssos_first"],
    minimize=cs_tssos_first_outer_dense["minimize_before_cs_tssos_first"],
    TS=cs_tssos_first_outer_dense["TS_before_cs_tssos_first"],
    merge=cs_tssos_first_outer_dense["merge_before_cs_tssos_first"],
    md=cs_tssos_first_outer_dense["md_before_cs_tssos_first"],
    QUIET=cs_tssos_first_outer_dense["QUIET_before_cs_tssos_first"],
    solver=cs_tssos_first_outer_dense["solver_before_cs_tssos_first"],
    tune=cs_tssos_first_outer_dense["tune_before_cs_tssos_first"],
    dualize=cs_tssos_first_outer_dense["dualize_before_cs_tssos_first"],
    solve=cs_tssos_first_outer_dense["solve_before_cs_tssos_first"],
    solution=cs_tssos_first_outer_dense["solution_before_cs_tssos_first"],
    Gram=cs_tssos_first_outer_dense["Gram_before_cs_tssos_first"],
    MomentOne=cs_tssos_first_outer_dense["MomentOne_before_cs_tssos_first"],
    Mommat=cs_tssos_first_outer_dense["Mommat_before_cs_tssos_first"],
    tol=cs_tssos_first_outer_dense["tol_before_cs_tssos_first"],
    cosmo_setting=cs_tssos_first_outer_dense["cosmo_setting_before_cs_tssos_first"],
    mosek_setting=cs_tssos_first_outer_dense["mosek_setting_before_cs_tssos_first"],
    normality=cs_tssos_first_outer_dense["normality_before_cs_tssos_first"],
    NormalSparse=cs_tssos_first_outer_dense["NormalSparse_before_cs_tssos_first"])


cs_tssos_first_inner_dense = deserialize("./example/data/nblockmix/cs_tssos_first_inner.jls")

cs_tssos_first_inner_dense["basis_if_initially_non_empty"]
cs_tssos_first_inner_dense["hbasis_if_initially_non_empty"]


cs_tssos_first_inner_dense["supp_before_resort"]
cs_tssos_first_inner_dense["coe_before_resort"]


cs_tssos_first_inner_dense["cliques_if_not_manually_specified"]
cs_tssos_first_inner_dense["cql_if_not_manually_specified"]
cs_tssos_first_inner_dense["cliquesize_if_not_manually_specified"]


cs_tssos_first_inner_dense["I_after_assign_constraint"]
cs_tssos_first_inner_dense["J_after_assign_constraint"]
cs_tssos_first_inner_dense["ncc_after_assign_constraint"]


cs_tssos_first_inner_dense["I_after_assign_constraint"]
cs_tssos_first_inner_dense["J_after_assign_constraint"]
cs_tssos_first_inner_dense["ncc_after_assign_constraint"]


cs_tssos_first_inner_dense["rlorder_if_uniform_degree"]


cs_tssos_first_inner_dense["cliques[i]_get_sbasis"]
cs_tssos_first_inner_dense["rlorder[i]_get_sbasis"]
cs_tssos_first_inner_dense["basis[i][1]_get_sbasis"]

cs_tssos_first_inner_dense["basis_if_initially_non_empty"]
cs_tssos_first_inner_dense["hbasis_if_initially_non_empty"]
	
cs_tssos_first_inner_dense["blocks_after_get_blocks"]
cs_tssos_first_inner_dense["eblocks_after_get_blocks"]
cs_tssos_first_inner_dense["cl_after_get_blocks"]
cs_tssos_first_inner_dense["blocksize_after_get_blocks"]


solvesdp_dense = deserialize("./example/data/nblockmix/solvesdp.jil")


cs_tssos_first_inner_dense["blocksize_before_sdpsolve"]

solvesdp_dense["tsupp_after_sadd"]


