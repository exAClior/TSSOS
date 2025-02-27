#import "@preview/algorithmic:0.1.0"
#import algorithmic: algorithm

#algorithm({
  import algorithmic: *
  Function("cs_tssos_first", args: ("pop", "x", "d", "keyword arguments"), {
		Cmt[obtain more efficient way to represent monomials] 

		Assign([ n, supp, coe ], FnI[polys_info][pop, x] )

		State[]

		Cmt[ Merge monomials with same support in objective]

    Assign([ `supp[1]`, `coe[1]` ], FnI[resort][`supp[1]`, `coe[1]`])

		State[]

		If(cond: [`clique` is not specified], {
			Assign([`cliques`], FnI[clique_decomp][n , m, numeq, dc, supp])
		})

		State[]

		Cmt[Group Constraints according to variables and cliques]

		Assign([I , J , ncc],FnI[assign_constraint][m, numeq, supp, cliques, cql])

		State[]

		Cmt[Determine max degree of monomials in subproblem for each clique]

		If(cond: [`d` is "min"],{
			Assign([`rlorder`], [max degree in I , J for each clique])	
		})

		Else({
			Assign([`rlorder`], [d for all cliques])	
		})

		State[]

		If(cond: [custom `basis` not specified],{

			Cmt[If Correlative Sparse, need to loop over all cliques]

			Assign([basis], FnI[get_sbasis][`cliques`, `rlorder`])	

			Cmt[ #text(red)[ Why is this sparse?]]

			Assign([hbasis],FnI[get_sbasis][])
		})

		If(cond: [Term Sparsity is required],{
			Assign([ksupp], [ supp + square of all terms in supp])

			Cmt[ Apply Term Sparsity by chordal extension, #text(red)[We can improve this]]
			Assign([blocks, eblocks, cl, blocksize],
			FnI[get_blocks][I, J , supp, cliques, cql, ...])
		})

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

