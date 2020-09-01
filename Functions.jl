# differential equations formula for finding densities of each pop at next time step
function fonction(u, p, t)
    x = u[1]
    y = u[2:end]
    regul = (1-(sum(u))/(p.K))
    dx = (p.bx*x + sum(p.ei.*y))*regul - p.ux*x-p.c*sum(p.βy.*y)*x
    dy = p.by .* y .* regul .- p.ui .* y .+ p.c .* p.βy .* x .* y
    return vcat(dx, dy)
end

function run_simulation()
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
end
