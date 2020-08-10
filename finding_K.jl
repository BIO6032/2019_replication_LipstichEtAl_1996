using DifferentialEquations
using Plots

n_parasites = 100

#Figure 2
c = 4.0;
ux = 0.2;
ux1 = fill(0.2, n_parasites) #le ux et uy 1000 à cause de la β # une autre façon :[0.2 for x in 1:1000]
uy = rand(200:1000, n_parasites)/1000 #doit être toujours >= à ux
u1 = (uy - ux1)
βy = 3 * u1 ./ (u1.+ 1) # augmente avec la mortalité u
bx = 1.0;
by = 0.1;
#ey = bx - by
ey = fill(0.9, n_parasites)    # constant for fig 2 part 1

Y = zeros(Float64, length(ey))
Y[1] = 0.0;
X0 = 10.0;


# SANDRINES TEST
debut = 0.0
duree = 1000.0
fin = debut + duree
N = zeros(Float64, (n_parasites+1, (n_parasites-1)*Int(duree)+1))
new_U = vcat(X0, Y)
parameters = (bx = bx, βy = βy, ey = ey, c = c, K = 100.0, ux = ux, by = by, uy = uy)


run_simulation()

Np = N'

plot(Np[:,2:end], c=:grey, lw=0.4, alpha=0.4)
plot!(Np[:,1], c=:black, lw=5, leg=false)
plot!(sum(Np[:,2:end]; dims=2))
