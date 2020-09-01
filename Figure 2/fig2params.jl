using DifferentialEquations
using Plots
import Random

n_parasites = 100;

ux = 0.2;
ux1 = fill(0.2, n_parasites); #le ux et ui 1000 à cause de la β # une autre façon :[0.2 for x in 1:1000]
Random.seed!(1234);
ui = rand(200:1000, n_parasites)/1000;
Random.seed!(1235);
r1 = rand(Float64, n_parasites)
Random.seed!(1236);
r2 = rand(Float64, n_parasites)
Random.seed!(1237);
r3 = rand(0.0:1000, n_parasites)/1000
bx = 1.0;
by = bx .* r1 .* (1 .- r1 .* r2);
# V0 = by*ux./(bx*ui_avg)
V0 = (by .* ux) ./ (bx .* ui);
α = 1 .- V0;
βy = r1 .-(α.* by)/bx
ei = bx .* (1 .- r3) .* (1 .- (r1 .* r2))

Y = zeros(Float64, length(ei));
Y[1] = 1.0;
X0 = 10.0;


debut = 0.0;
duree = 1000.0;
fin = debut + duree;
N = zeros(Float64, (n_parasites+1, (n_parasites-1)*Int(duree)+1));
new_U = vcat(X0, Y);
parameters = (bx = bx, βy = βy, ei = ei, c = c, K = 100.0, ux = ux, by = by, ui = ui);
