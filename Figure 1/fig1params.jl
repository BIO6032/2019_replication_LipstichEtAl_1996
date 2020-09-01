using DifferentialEquations
using Plots
import Random


n_parasites = 100;

#parameters
c = 4.0;
ux = 0.2;
ux1 = fill(0.2, n_parasites); #le ux et ui 1000 à cause de la β # une autre façon :[0.2 for x in 1:1000]
Random.seed!(1234);
ui = rand(200:1000, n_parasites)/1000; #always >= ux
u1 = (ui - ux1);
βy = 3 * u1 ./ (u1.+ 1) # augmente avec la mortalité u
bx = 1.0;
ei = fill(0.9, n_parasites)   # constant for fig 2 part 1

debut = 0.0;
duree = 1000.0;
fin = debut + duree;
N = zeros(Float64, (n_parasites+1, (n_parasites-1)*Int(duree)+1));
parameters = (bx = bx, βy = βy, ei = ei, c = c, K = 99, ux = ux, by = by, ui = ui);

Y = zeros(Float64, length(ei));
X0 = 10.0;
new_U = vcat(X0, Y);
Y[1] = 1.0;
