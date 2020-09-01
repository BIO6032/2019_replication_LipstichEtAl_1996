bi = 1.0;
include("fig1params.jl")


include("../Functions.jl")
run_simulation()


############ Figure 1b ############

Np = N';

lbls = ["" for i = 1:1:n_parasites];
lbls2 = vcat("Infected", lbls);
lbls2 = hcat(lbls2...);

plot(Np[:,2:end], c=:blue, lw=0.4, alpha=0.4, title = "Number of infected and uninfected hosts",
    xlabel = "Time", ylabel = "Number of individuals", label = lbls2, ylims =(0,100))
plot!(Np[:,1], c=:black, lw=0.4, label = "Uninfected")
#plot!(sum(Np[:,2:end]; dims=2), label = "total parasites")

# png("Figure 1/graph_1b.png")


########### Figure 1h #############

#strains that survive at each time step for ui
survival = (Np.>0.0)[:,2:end];
survived_ui = survival.*ui';
avg_survived = sum(survived_ui; dims=2)./sum(survival; dims=2);

avg_w_survived = sum(Np[:,2:end].*ui'; dims=2)./sum(Np[:,2:end]; dims=2);
ui_avg = avg_w_survived;

# plot(avg_survived, title = "Average ui in the population",xlabel = "Time", ylabel = "Mean mortality (ui)", leg = false)
plot(avg_w_survived, c=:black,title = "Average mortality in the population",
    xlabel = "Time", ylabel = "Mean mortality (ui)", leg = false, ylims =(0,1))

# png("Figure 1/graph_1h.png")


############ Figure 2f ##############

#strains that survive for βy
survived_βy = survival.*βy';
avg_survived_βy = sum(survived_βy; dims=2)./sum(survival; dims=2);
βi_avg = avg_survived_βy;

#avec moins de bruit (weighted)
avg_w_survived_βy = sum(Np[:,2:end].*βy'; dims=2)./sum(Np[:,2:end]; dims=2);
βy_avg = avg_w_survived_βy;

#calculating R0
k=1;
H0 = c*βi_avg./ui_avg.*k.*(1-ux/bx);

H0_w = c*βy_avg./ui_avg.*k.*(1-ux/bx);

V0 = bi*ux./(bx*ui_avg);

plot(V0,c=:black, title = "Average vertical transmission in the population", xlabel = "Time", ylabel = "Mean V0", leg = false, ylims=(0,1))

# png("Figure 1/graph_1f")


######## Figure 1d ##########

R0 = H0 + V0;

R0_w = H0_w + V0;

# plot(R0, title = "Average R0 in the population", xlabel = "Time", ylabel = "Mean R0", leg = false, ylims=(0,))
plot(R0_w,c=:black, title = "Average R0 in the population \nwith R0 = V0 +H0", xlabel = "Time", ylabel = "Mean R0", leg = false, ylims=(0,7))

# png("Figure 1/graph_1d")



########## Figure 1j ###########

#calculating evenness
function pielou(n)
    np = filter(x -> x > eps(), n)
    p = np./sum(np)
    ev = length(p) == 1 ? 0.0 : -sum(p.*log.(p))*(1/log(length(p)))
    return ev
end;
ev = mapslices(pielou, Np[:,2:end]; dims=2);
plot(ev, c=:black, title = "Evenness", xlabel = "Time", ylabel = "Relative abundance (log)", leg = false, ylims = (0,1))

# png("Figure 1/graph_1j.png")
