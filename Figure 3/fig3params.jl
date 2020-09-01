using DifferentialEquations
using Plots
using Distributions
import Random

n_parasites = 100;

ux = 0.2;
ux1 = fill(0.2, n_parasites); #le ux et uy 1000 à cause de la β # une autre façon :[0.2 for x in 1:1000]
Random.seed!(1234);
uy = rand(200:1000, n_parasites)/1000;
Random.seed!(1235);
r1 = rand(Float64, n_parasites);
Random.seed!(1236);
r2 = rand(Float64, n_parasites);
Random.seed!(1237);
r3 = rand()*r1
bx = 1.0;
by = bx .* r3 .* (1 .- r1 .* r2);
V0 = (by .* ux) ./ (bx .* uy);
α = 1 .- V0;
βy = r1 .-(α.* by)/bx;
ey = bx .* (1 .- r3) .* (1 .- (r1 .* r2));

Y = zeros(Float64, length(ey));
Y[1] = 1.0;
X0 = 10.0;

debut = 0.0;
duree = 1000.0;
fin = debut + duree;
N = zeros(Float64, (n_parasites+1, (n_parasites-1)*Int(duree)+1));
new_U = vcat(X0, Y);
parameters = (bx = bx, βy = βy, ey = ey, c = c, K = 100.0, ux = ux, by = by, uy = uy);
