#calculating evenness
include("fig3b.jl")

function pielou(n)
    np = filter(x -> x > eps(), n)
    p = np./sum(np)
    ev = length(p) == 1 ? 0.0 : -sum(p.*log.(p))*(1/log(length(p)))
    return ev
end
ev = mapslices(pielou, Np[:,2:end]; dims=2);
plot(ev,c=:black, title = "Evenness", xlabel = "Time", ylabel = "Relative abundance (log)", leg = false,ylims = (0,1))

# png("Figure 3/graph_3j.png")
