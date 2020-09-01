using DifferentialEquations
using Plots
using Distributions

n_parasites = 100

c = 4.0;
ux = 0.2;
ux1 = fill(0.2, n_parasites); #le ux et uy 1000 à cause de la β # une autre façon :[0.2 for x in 1:1000]
Random.seed!(1234);
uy = rand(200:1000, n_parasites)/1000;
Random.seed!(1235);
r1 = rand(Float64, n_parasites);
Random.seed!(1236);
r2 = rand(Float64, n_parasites);
Random.seed!(1237);
r3 = rand()*r1;
bx = 1.0;
by = bx .* r1 .* (1 .- r1 .* r2);
V0 = (by .* ux) ./ (bx .* uy);
α = 1 .- V0;
βy = r1 .-(α.* by)/bx;
ey = bx .* (1 .- r3) .* (1 .- (r1 .* r2));

Y = zeros(Float64, length(ey));
Y[1] = 1.0;
X0 = 10.0;

# differential equations formula for finding densities of each pop at next time step
function fonction(u, p, t)
    x = u[1]
    y = u[2:end]
    regul = (1-(sum(u))/(p.K))
    dx = (p.bx*x + sum(p.ey.*y))*regul - p.ux*x-p.c*sum(p.βy.*y)*x
    dy = p.by .* y .* regul .- p.uy .* y .+ p.c .* p.βy .* x .* y
    return vcat(dx, dy)
end

debut = 0.0;
duree = 1000.0;
fin = debut + duree;
N = zeros(Float64, (n_parasites+1, (n_parasites-1)*Int(duree)+1));
new_U = vcat(X0, Y);
parameters = (bx = bx, βy = βy, ey = ey, c = c, K = 80.0, ux = ux, by = by, uy = uy);

# each strain introduction (1000x)
@progress "Simulation" for i in 2:length(Y)
    # initial conditions
    prob = ODEProblem(fonction, new_U, (debut,fin), parameters)
    solution = solve(prob, saveat=debut:1.0:fin)
    for t in eachindex(solution.t)
        pop = solution.u[t]
        for i in 1:100
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

#strains that survive at each time step
survival = (Np.>0.0)[:,2:end];

#mean β
survived_βy = survival.*βy';
avg_survived_βy = sum(survived_βy; dims=2)./sum(survival; dims=2);
βy_avg = avg_survived_βy;

#weighted β
avg_w_survived_βy = sum(Np[:,2:end].*βy'; dims=2)./sum(Np[:,2:end]; dims=2);
βi_avg = avg_w_survived_βy;

#mean virulence
Vir = (1 .- (by .+ ey).*ux./(bx.*uy));

survived_Vir = survival.*Vir';
avg_survived_Vir = sum(survived_Vir; dims=2)./sum(survival; dims=2);
Vir_avg = avg_survived_Vir;

#weighted virulence
avg_w_survived_Vir = sum(Np[:,2:end].*Vir'; dims=2)./sum(Np[:,2:end]; dims=2);
Vir_avg_w = avg_w_survived_Vir;

#plot(βy_avg, title = "Mean virulence and beta", label = "beta")
#plot!(Vir_avg, label = "virulence")
plot(βi_avg,c=:black, title = "Virulence and beta", label = "Beta", xlabel = "Time", ylabel = "Mean virulence & \n Mean Beta", ylims = (0,1))
plot!(Vir_avg_w,c=:blue, label = "Virulence")

# png("Figure 4/graph_4h.png")
