#import "@preview/algorithmic:0.1.0"
#import algorithmic: algorithm

#algorithm({
  import algorithmic: *
  Function("cs_tssos_first", args: ("pop", "x", "d", "keyword arguments"), {
		Cmt[obtain more efficient way to represent monomials] 

		Assign([ n, supp, coe ], FnI[ polys_info][pop, x] )

		Cmt[ Merge monomials with same support in objective]

    Assign([supp[1], coe[1]],FnI[resort][supp[1], coe[1]])) 

		Return[*null*]
  })
})



#algorithm({
  import algorithmic: *
  Function("cs_tssos_first", args: ("A", "n", "v"), {
    Cmt[Initialize the search range]
    Assign[$l$][$1$]
    Assign[$r$][$n$]
    State[]
    While(cond: $l <= r$, {
      Assign([mid], FnI[floor][$(l + r)/2$])
      If(cond: $A ["mid"] < v$, {
        Assign[$l$][$m + 1$]
      })
      ElsIf(cond: [$A ["mid"] > v$], {
        Assign[$r$][$m - 1$]
      })
      Else({
        Return[$m$]
      })
    })
    Return[*null*]
  })
})

#algorithm({
  import algorithmic: *
  Function("tssos_first", args: ("A", "n", "v"), {
    Cmt[Initialize the search range]
    Assign[$l$][$1$]
    Assign[$r$][$n$]
    State[]
    While(cond: $l <= r$, {
      Assign([mid], FnI[floor][$(l + r)/2$])
      If(cond: $A ["mid"] < v$, {
        Assign[$l$][$m + 1$]
      })
      ElsIf(cond: [$A ["mid"] > v$], {
        Assign[$r$][$m - 1$]
      })
      Else({
        Return[$m$]
      })
    })
    Return[*null*]
  })
})

matrixsos.jl LinearPMI_first