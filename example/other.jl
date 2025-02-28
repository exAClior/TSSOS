using TSSOS, DynamicPolynomials

@polyvar x[1:3]
f = 1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2]*x[3] + x[2]
g = 1 - x[1]^2 - 2*x[2]^2
h = x[2]^2 + x[3]^2 - 1
pop = [f, g, h]
d = 2 # set the relaxation order
opt,sol,data = tssos_first(pop, x, d, numeq=1, TS="MD") # compute the first TS step of the TSSOS hierarchy