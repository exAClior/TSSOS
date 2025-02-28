#import "@preview/algorithmic:0.1.0"
#import algorithmic: algorithm

#text(red)[ We did not consider Term Sparsity in this algorithm]

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

		Cmt[#text(red)[When do we have non-empty `constraint_not_belong_to_single_clique`?]]
		Assign([`clique_eq_cons` , `clique_ineq_cons`, `constraint_not_belong_to_single_clique`],
           FnI[assign_constraint][`eq_cons`, `ineq_cons`, `cliques`])

		State[]

		Cmt[If Correlative Sparse, need to loop over all cliques]

		Cmt[Get Standard Monomial Basis upto some order] 

		Assign([`obj_basis`], FnI[get_sbasis][`clique_eq_cons`, $d$])

		State[]

		Assign([`halfdegree_of_each_constraint`],FnI[halfdegree][`clique_ineq_cons`])

		Assign([`cons_basis`],FnI[get_sbasis][`clique_ineq_cons`, $d$ - `halfdegree_of_each_constraint`])

		State[]

		Assign([`opt`],FnI[solvesdp][$n$,`obj`, `eq_cons`, `ineq_cons`, `obj_basis`, `cons_basis`])

		State[]

		Return[`opt`]
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

$1_(alpha,bold(0)) = cases(
  1 "if" alpha = bold(0),
  0 "otherwise"
)$




#algorithm({
  import algorithmic: *
  Function("solvesdp", args: ([$n$], [`obj`], [`eq_cons`], [`ineq_cons`], [`obj_basis`], [`cons_basis`]), {

		Assign([`tsupp`], [#FnI[get_support][`obj_basis`] $union$ ($union_(sigma in #(`cons_basis`))$ #FnI[get_support][$sigma$]) ])
		Assign([`all_cons`],[`{1}` $union$ `ineq_cons` $union$ `eq_cons`]) 

		Assign([`b`],FnI[init_variable][1])

		State[]

		Assign([`G_matrices`], [( )]) // Initialize an empty collection of matrices
		Assign([`j`], [0])
		
		For(cond:[`cur_basis` in (`obj_basis` $union$ `cons_basis`)],{
		    Assign([`G_matrices[j]`], FnI[init_psd_matrix_variable][length(`cur_basis`)])
		    Assign([`j`], [`j` + 1])
		})

		State[]
		
		Assign([`SDP_constraints`], FnI[init_sdp_constraints][])

		For(cond:[`alpha` in `tsupp`],{
			Cmt[Get coefficient of monomial in objective]
			Assign([`f_alpha`], FnI[get_coefficient][`obj`, `alpha`])

			Assign([`include_b`], FnI[is_equal][`alpha`, $bold(0)$])	

			Assign([`constraint_sum`], [0]) // Initialize sum of traces
			Assign([`j`], [0])
			For(cond:[`cur_basis` in (`obj_basis` $union$ `cons_basis`)],{
			    Assign([`c_alpha_j`], FnI[compute_c_alpha_j][`alpha`, $#(`all_cons`)_#(`j`)$, `cur_basis`])
			    Assign([`trace_j`], FnI[trace][`c_alpha_j` $dot$ `G_matrices[j]`])
			    Assign([`constraint_sum`], [`constraint_sum` + `trace_j`])
			    Assign([`j`], [`j` + 1])
			})
			
			Assign([`SDP_cons_alpha`], FnI[set_constraint_equal][`f_alpha` - `b` $dot$ `include_b`, `constraint_sum`])
			Assign([`SDP_constraints`], [`SDP_constraints` $union$ `SDP_cons_alpha`])
		})
		
		Assign([`opt`],FnI[optimize][`G_matrices`, `b`, `SDP_constraints`])

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
  Function("compute_c_alpha_j", args: ([`alpha`],[`constraint_poly`] ,[`cur_basis`]), {
		Assign([`c_alpha_j`], FnI[init_matrix][length(`cur_basis`), length(`cur_basis`)]) // Initialize matrix
		
		For(cond:[ i in 1:length(`cur_basis`)],{
			For(cond:[ k in 1:length(`cur_basis`)],{
				Assign([`cur_polynomial`], [ $#(`constraint_poly`) dot #(`cur_basis`)_i dot #(`cur_basis`)_k$])
				If(cond: [`alpha` in `cur_polynomial` ],{
					Assign([$#(`c_alpha_j`)_(i,k)$], FnI[get_coefficient][`cur_polynomial`, `alpha`])
				})
				Else({
					Assign([$#(`c_alpha_j`)_(i,k)$], [0])
				})
			})
		})
		
		Return[`c_alpha_j`] // Return the matrix
  })
})

#algorithm({
  import algorithmic: *
  Function("halfdegree", args: ([`constraints`],), {
    Assign([`half_degrees`], [$emptyset$])
    For(cond:[`constraint` in `constraints`],{
        Assign([`max_degree`], FnI[get_max_degree][`constraint`])
        Assign([`half_degree`], FnI[ceil][`max_degree` / 2])
        Assign([`half_degrees`], [`half_degrees` $union$ `half_degree`])
    })
    Return[`half_degrees`]
  })
})

#algorithm({
  import algorithmic: *
  Function("init_sdp_constraints", args: (), {
    Return[{}] // Empty set of constraints
  })
})

// #algorithm({
//   import algorithmic: *
//   Function("approx_sol", args: ([`opt`], [`moment`], [`n`], [`clique_eq_cons`], [`clique_ineq_cons`], [`obj`], [`eq_cons`], [`ineq_cons`]), {

//   })
// })


#algorithm({
  import algorithmic: *
  Function("cs_nctssos_first", args: ("obj", "eq_cons", "ineq_cons", $n$, $d$), {

		Cmt[ Canonicalize the support of objective ]
		Assign([ `obj_supp`, `obj_coe` ], FnI[canonicalize][`obj`]) // Fixed arguments
		
		State[]
		
		Assign([`cliques`], FnI[clique_decomp][$n$, `obj`, `eq_cons`, `ineq_cons`])
		
		State[]
		
		Cmt[Get basis for objective and constraints]
		Assign([`all_obj_basis`], [{}])
		Assign([`all_cons_basis`], [{}])
		
		For([`clique` in `cliques`],{
		    Assign([`clique_size`], FnI[length][`clique`]) // Define clique_size
		    
		    Assign([`obj_basis_for_clique`], FnI[get_ncbasis][`clique_size`, `d`, `clique`])
		    Assign([`all_obj_basis`], [`all_obj_basis` $union$ `obj_basis_for_clique`])
		    
		    Assign([`cons_basis_for_clique`], FnI[get_ncbasis][`clique_size`, `d`, `clique`])
		    Assign([`all_cons_basis`], [`all_cons_basis` $union$ `cons_basis_for_clique`])
		    
		    If(cond: [FnI[has_commutative_vars][`clique`]],{
		        Assign([`obj_basis_for_clique`], FnI[simplify][`obj_basis_for_clique`])
		        Assign([`cons_basis_for_clique`], FnI[simplify][`cons_basis_for_clique`])
		    })
		})
		
		Cmt[Reduce the support by considering commutative portion of variables]
		Assign([`tsupp`], FnI[get_support][`all_obj_basis` $union$ `all_cons_basis`])
		Assign([`tsupp`], FnI[reduce_support][`tsupp`, `obj`])
		
		Cmt[ Apply Term Sparsity by chordal extension]
		Assign([`blocks`, `eblocks`, `cl`, `blocksize`],
		    FnI[get_blocks][`tsupp`, `cliques`])
		
		Assign([`opt`, `moment`, `GramMat`],FnI[solvesdp][$n$, `obj`, `eq_cons`, `ineq_cons`, 
		                                         `all_obj_basis`, `all_cons_basis`])
		
		Cmt[Construct the data structure for higher order improvement]
		Assign([`data`], FnI[mcpop_data][$n$, `blocks`, `eblocks`, `cl`, `blocksize`])
		
		Return[`opt`, `moment`, `data`]
  })
})

