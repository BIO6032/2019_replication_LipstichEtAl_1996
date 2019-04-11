using DifferentialEquations
using Plots

#Figure 2
c = 4.0
ux = 0.2
bx = 1.0
by = 0.1
uy = rand(Float64)

# struct hote
#     ux :: Float64
#     bx :: Float64
# end

# struct parasites
#     uy :: Float64
#     by :: Float64
#     #ey :: Float64
#     #βy :: Float64
# end

# function new_para()
#     uy = rand(Float64)
#     by = 0.1
#     new_Y = parasites(uy, by)
#     return new_Y
# end

#newY = new_para()

ey = fill(0.9, 100)    # constant for fig 2 part 1
βy = fill(0.0, 100)

Y = zeros(Float64, length(ey))
Y[1] = 10.0
X0 = 1.0

function fu(u, p, t)
    x = u[1]
    y = u[2:end]
    regul = (1-(sum(u))/(p.K))
    dx = (p.bx*x + sum(p.ey.*y))*regul - p.ux*x-p.c*sum(p.βy.*y)*x
    dy = p.by .* y .* regul .- p.uy .* y .+ p.c .* p.βy .* x .* y
    return vcat(dx, dy)
end

parameters = (bx = bx, βy = βy, ey = ey, c = c, K = 20.0, ux = ux, by = by, uy = uy)

prob = ODEProblem(fu, vcat(X0, Y), (0.0, 100.0), parameters)
solution = solve(prob)
plot(solution, leg=false)

#new_u = solution[end]

# new_y = findfirst(x -> x == 0.0, new_u)
# new_u[new_y] = 1
#
# prob = ODEProblem(fu, new_u, (100.0, 200.0), parameters)
# solution = solve(prob)
#
# plot(solution, leg=false)
