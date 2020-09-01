c = 0.5;
include("fig2params.jl");


include("../Functions.jl")
run_simulation()


########## Figure 2a ##########
Np = N';

lbls = ["" for i = 1:1:n_parasites];
lbls2 = vcat("Infected", lbls);
lbls2 = hcat(lbls2...);

plot(Np[:,2:end], c=:blue, lw=0.4, alpha=0.4, title = "Number of infected and uninfected hosts",
    xlabel = "Time", ylabel = "Number of individuals", label = lbls2, ylims=(0,100))
plot!(Np[:,1], c=:black, lw=0.4, label = "Uninfected")
# plot!(sum(Np[:,2:end]; dims=2), label = "Total # parasites")

# png("Figure 2/graph_2a.png")


######### Figure 2e ##########

#strains that survive at each time step for uy
survival = (Np.>0.0)[:,2:end];
survived_ui = survival.*uy';
avg_survived = sum(survived_ui; dims=2)./sum(survival; dims=2);
avg_w_survived = sum(Np[:,2:end].*uy'; dims=2)./sum(Np[:,2:end]; dims=2);
uy_avg = avg_w_survived;

#strains that survive for βy
survived_βy = survival.*βy';
avg_survived_βy = sum(survived_βy; dims=2)./sum(survival; dims=2);
βi_avg = avg_survived_βy;

#weighted β
avg_w_survived_βy = sum(Np[:,2:end].*βy'; dims=2)./sum(Np[:,2:end]; dims=2);
βy_avg = avg_w_survived_βy;

survived_by = survival.*by';
avg_by_survived = sum(survived_by; dims=2)./sum(survival; dims=2);;
avg_w_by = sum(Np[:,2:end].*by'; dims=2)./sum(Np[:,2:end]; dims=2)
by_avg = avg_w_by;

V0= by_avg.*ux./(bx*uy_avg);
plot(V0,c=:black, title = "Vertical cases", xlabel = "Time", ylabel = "Mean V0", leg = false, ylims = (0.0,1.0))
# png("Figure 2/graph_2e")


########## Figure 2c ############

#calculating R0
k=1;
H0 = c*βi_avg./uy_avg.*k.*(1-ux/bx);

H0_w = c*βy_avg./uy_avg.*k.*(1-ux/bx);

V0_w = by_avg.*ux./(bx*uy_avg);

R0 = H0 + V0_w;

R0_w = H0_w + V0_w;

# plot(R0, title = "Average R0 in the population", xlabel = "Time", ylabel = "Mean R0", leg = false)
plot(R0_w,c=:black, title = "Average R0 in the population \nwith R0 = V0 +H0", xlabel = "Time", ylabel = "Mean R0", leg = false, ylims=(0,7))

# png("Figure 2/graph_2c")



######### Figure 2g ###########

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

# png("Figure 2/graph_2g.png")


########## Figure 2i ##########

function pielou(n)
    np = filter(x -> x > eps(), n)
    p = np./sum(np)
    ev = length(p) == 1 ? 0.0 : -sum(p.*log.(p))*(1/log(length(p)))
    return ev
end
ev = mapslices(pielou, Np[:,2:end]; dims=2);
plot(ev,c=:black, title = "Evenness", xlabel = "Time", ylabel = "Relative abundance (log)", leg = false,ylims = (0,1))

# png("Figure 2/graph_2i.png")
