using DifferentialEquations
using Plots
import Random


n_parasites = 100;

#parameters
c = 4.0;
ux = 0.2;
ux1 = fill(0.2, n_parasites); #le ux et uy 1000 à cause de la β # une autre façon :[0.2 for x in 1:1000]
Random.seed!(1234);
uy = rand(200:1000, n_parasites)/1000; #always >= ux
u1 = (uy - ux1);
βy = 3 * u1 ./ (u1.+ 1); # augmente avec la mortalité u
bx = 1.0;
by = 0.1;
ey = fill(0.9, n_parasites);    # constant for fig 2 part 1
debut = 0.0;
duree = 1000.0;
fin = debut + duree;
N = zeros(Float64, (n_parasites+1, (n_parasites-1)*Int(duree)+1));
parameters = (bx = bx, βy = βy, ey = ey, c = c, K = 80.0, ux = ux, by = by, uy = uy);

Y = zeros(Float64, length(ey));
X0 = 10.0;
new_U = vcat(X0, Y);
Y[1] = 1.0;


# differential equations formula for finding densities of each pop at next time step
function fonction(u, p, t)
    x = u[1]
    y = u[2:end]
    regul = (1-(sum(u))/(p.K))
    dx = (p.bx*x + sum(p.ey.*y))*regul - p.ux*x-p.c*sum(p.βy.*y)*x
    dy = p.by .* y .* regul .- p.uy .* y .+ p.c .* p.βy .* x .* y
    return vcat(dx, dy)
end

# each strain introduction (1000x)
@progress "Simulation" for i in 2:length(Y)
    # initial conditions
    prob = ODEProblem(fonction, new_U, (debut,fin), parameters)
    solution = solve(prob, saveat=debut:1.0:fin)
    for t in eachindex(solution.t)
        pop = solution.u[t]
        for i in 1:n_parasites
            if (pop.<0)[i]
                pop[i] = 0
            end
        end
        N[:,Int(solution.t[t]+1)] = pop
    end

    # set conditions & new parasite for next loop
    global new_U = solution[end]
    new_y = findfirst(x -> x == 0.0, new_U)
    new_U[new_y] = 1.0

    # set limits for next loop
    global debut = fin
    #global fin = Float64((i+1).* 1000.0)
    global fin = debut + duree
end

Np = N';

lbls = ["" for i = 1:1:n_parasites];
lbls2 = vcat("Parasites", lbls);
lbls2 = hcat(lbls2...);

plot(Np[:,2:end], c=:blue, lw=0.4, alpha=0.4, title = "Number of infected and uninfected hosts",
    xlabel = "Time", ylabel = "Number of individuals", label = lbls2, ylims =(0,70))
plot!(Np[:,1], c=:black, lw=0.4, label = "Hosts")
#plot!(sum(Np[:,2:end]; dims=2), label = "total parasites")

# png("Figure 2/graph_2a.png")
