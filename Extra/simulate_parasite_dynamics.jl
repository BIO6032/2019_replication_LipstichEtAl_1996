using Plots
include("Fig1_strain_parameters.jl")

#using fonctions defined in Fig1_strain_parameters
# making a new parasite & a new host
p = new_Y()
h = new_X()
#model
N = parasites(h, p, 0.5, 0.5, timesteps = 2000, iter = 5000)
#graph
plot(N, labels = ["X", "Y"], title = "Dynamics of populations X and Y")

png("img/one_parasite.png")
