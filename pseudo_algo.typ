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

		Cmt[#text(red)[When do we have non-empty `contraint_not_belong_to_single_clique`?]]
		Assign([`clique_eq_cons` , `clique_ineq_cons`, `constraint_not_belong_to_single_clique`],FnI[assign_constraint][`eq_cons`, `ineq_cons`, `cliques`])

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

$
	op("sup",limits:#true)_(bold(G)_j,b) b	\
	s.t. f_(alpha)	- b 1_(alpha,bold(0)) = sum_(j = 0)^(m) angle.l bold(C)_(alpha)^(j), bold(G)_j angle.r, alpha in bb(N)_(2r)^(n) \
	 bold(G)_j in bb(S)_(t_j)^(+) , j = 0, ..., m
$

where $t_j$ is the moment matrix size of $j$th constraint or objective.

$1_(alpha,bold(0))$ is $0$ unless $alpha = bold(0)$.




#algorithm({
  import algorithmic: *
  Function("solvesdp", args: ([`n`], [`obj`], [`eq_cons`], [`ineq_cons`], [`obj_basis`], [`cons_basis`]), {

		Cmt[Get monomials appeared in moment matrices]
		Assign([`tsupp`], [#FnI[get_support][`obj_basis`] $union$ ($union_(sigma in #(`cons_basis`))$ #FnI[get_support][$sigma$]) ])
		Assign([`all_cons`],[`{1}` $union$ `ineq_cons` $union$ `eq_cons`]) 

		Assign([`b`],FnI[init_variable][1])
		For(cond:[(`j`, `cur_basis`) in enumerate(`cons_basis` $union$ `obj_basis`)],{
			Assign([$G_j$], FnI[init_variable][length(`cur_basis`)])
		})

		For(cond:[`alpha` in `tsupp`],{
			Cmt[Get coefficient of monomial in objective]
			Assign([`f_alpha`], FnI[get_coefficient][`obj`, `alpha`])

			Assign([`include_b`], FnI[is_equal][`alpha`, $bold(0)$])	

			For(cond:[(`j`, `cur_basis`) in enumerate(`obj_basis` $union$ `cons_basis`)],{
				Assign([`c_alpha_j`], FnI[compute_c_alpha_j][`alpha`, $#(`all_cons`)_j$ ,`cur_basis`])
			})

			Assign([`j_trace`], FnI[trace][`c_alpha_j` $dot$ `G_j`])
			Assign([`SDP_cons_j`], FnI[set_constraint_equal][`f_alpha` - b $dot$ `include_b`, `j_trace`])

		})
	Assign([`opt`],FnI[optimize][$G_j$, `b`, `SDP_cons_j`])

	Return([`opt`])
  })
})



#algorithm({
  import algorithmic: *
  Function("get_support", args: ([`basis`],), {
		Assign([`supports`], [$emptyset$])
		For(cond:[i in `basis`],{
			For(cond:[j in `basis`],{
				Assign([`supports`], [`supports` $union$ (i $dot$ j)])
			} )
		})
		Return[`supports`]
  })
})

#algorithm({
  import algorithmic: *
  Function("compute_c_alpha_j", args: ([`alpha`],[`constraint`] ,[`cur_basis`]), {
		For(cond:[ i in 1:length(`cur_basis`)],{
			For(cond:[ k in 1:length(`cur_basis`)],{
				Assign([`cur_term`], [ $#(`constraint`) dot #(`cur_basis`)_i dot #(`cur_basis`)_k$])
				If(cond: [`alpha` in `cur_term` ],{
					Assign([$#(`c_alpha_j`)_(i,k)$], FnI[get_coefficient][`alpha`, `cur_term`])
				})
				Else({
					Assign([$#(`c_alpha_j`)_(i,k)$], [0])
				})
			})
		})
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

