#import "@preview/algorithmic:0.1.0"
#import algorithmic: algorithm

#text(red)[ We did not consider Term Sparsity in this algorithm]

- pop: list of polynomials containing objective and constraints
- n: number of variables
- supp: support of the monomials
- coe: coefficients of the monomials
- m: number of constraints
- numeq: number of equality constraints
- dc: degree of the constraints and objectives 

#algorithm({
  import algorithmic: *
  Function("cs_tssos_first", args: ([`obj`],[`eq_cons`],[`ineq_cons`], [$n$], [$d$]), {

		Cmt[ Merge monomials with same support in objective]

    Assign([ `obj`], FnI[resort][`obj`])

		State[]

		Cmt[Automatically generate cliques]
		Assign([`cliques`], FnI[clique_decomp][$n$ , `obj`, `eq_cons`, `ineq_cons`])

		State[]

		Cmt[Group Constraints according to variables and cliques]

		Assign([`clique_eq_cons` , `clique_ineq_cons`, `clique_with_no_constraint`],FnI[assign_constraint][`eq_cons`, `ineq_cons`, `cliques`])

		State[]

		Cmt[If Correlative Sparse, need to loop over all cliques]

		Cmt[Get Standard Monomial Basis upto some order] 

		Assign([`obj_basis`], FnI[get_sbasis][`cliques_eq_cons`, $d$])	

		State[]

		Assign([`halfdegree_of_each_constraint`],FnI[halfdegree][[`cliques_ineq_cons`]])

		Assign([`cons_basis`],FnI[get_sbasis][`cliques_ineq_cons`, $d$ - `halfdegree_of_each_constraint`])

		State[]

		Assign([`opt`, `moment`],FnI[solvesdp][$n$,`obj`, `eq_cons`, `ineq_cons`, `obj_basis`, `cons_basis`])

		State[]

		Cmt[Why is this step necessary?]

		Assign([sol], FnI[approx_sol][`opt`, `moment`, $n$, `clique_eq_cons`, `clique_ineq_cons`, `obj`, `eq_cons`, `ineq_cons`])

		Return[opt, sol]
  })
})

#algorithm({
  import algorithmic: *
  Function("clique_decomp", args: ([$n$], [`obj`], [`eq_cons`], [`ineq_cons`]), {

		Assign([`G`],FnI[get_variable_dependency_graph][$n$,`obj`, `eq_cons`, `ineq_cons`])

		Assign([`G`], FnI[chordal_extention][`G`])

		Assign([`cliques`], FnI[get_cliques][`G`])

		Return[`cliques`]
  })
})


#algorithm({
  import algorithmic: *
  Function("solvesdp", args: ([`n`], [`obj`], [`eq_cons`], [`ineq_cons`], [`obj_basis`], [`cons_basis`]), {

		Assign([`tsupp`], FnI[get_monomials_in_moment_mat][`obj_basis`, `cons_basis`])
		

  })
})

#algorithm({
  import algorithmic: *
  Function("approx_sol", args: ([`opt`], [`moment`], [`n`], [`clique_eq_cons`], [`clique_ineq_cons`], [`obj`], [`eq_cons`], [`ineq_cons`]), {

  })
})


#algorithm({
  import algorithmic: *
  Function("cs_nctssos_first", args: ("pop", "x", "d", "keyword arguments"), {
		Cmt[obtain more efficient way to represent monomials] 

		Assign([n, supp, coe ], FnI[polys_info][pop, x] )

		State[]

		Cmt[ Merge terms with canonically same support in objective]

		If(cond: [`obj` is "trace"],{
			Assign([ `obj_supp`, `obj_coe` ], FnI[cyclic_canon][`obj_supp`, `obj_coe`])
		})
		ElsIf(cond: [`obj` is "eigen"],{
			Assign([ `obj_supp`, `obj_coe` ], FnI[sym_canon][`obj_supp`, `obj_coe`])
		})

		State[]

		Cmt[Decompose into cliques, we can come in]
		If(cond: [Require Correlative Sparsity], {
			Assign([`cliques`], FnI[clique_decomp][n , m, numeq, dc, supp])
		})

		State[]

		Cmt[Group Constraints according to variables and cliques]

		Assign([J , ncc],FnI[assign_constraint][m, supp, cliques, cql])

		State[]


		Cmt[Get basis for objective and constraints]

		For([clique in `cliques`],{

			Assign([objective basis], FnI[get_ncbasis][`cliquesize`, `d`, `clique`])	

			Assign([constraint basis], FnI[get_ncbasis][`cliquesize`, `d`, `clique`])	

			If(cond: [variables contains commutative ones],{
				Assign([basis], FnI[simplify][`objective basis`, `constraint basis`])
			})
		})


		Cmt[Reduce the support by considering commutative portion of variables]
		Assign([tsupp], FnI[reduce!.][tsupp, obj=obj, par])

		Cmt[ Apply Term Sparsity by chordal extension, #text(red)[We can improve this]]
		Assign([blocks, eblocks, cl, blocksize],
			FnI[get_blocks][I, J , supp, cliques, cql, ...])

		Assign([opt, ksupp, moment, GramMat ...],FnI[solvesdp][n,m,...])

		Cmt[Construct the data structure for higher order improvement]
		Assign([data], FnI[mcpop_data][n, nb, m, numeq, ...])

		Cmt[Why is this step necessary?]
		If(cond: [`solution` is true],{
			Assign([sol, gap, data.flag], FnI[approx_sol][n,m,...])
		})

		Return[opt, sol, data]
  })
})

