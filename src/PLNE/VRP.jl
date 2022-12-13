using JuMP
using CPLEX
using CPUTime
using Graphs
import HiGHS

include("../utils.jl")

function VRP(prp, t, optim="CPLEX")
    if optim == "HiGHS"
        LP = Model(HiGHS.Optimizer)
    elseif optim == "CPLEX"
        LP = Model(CPLEX.Optimizer)
    else
        return -1
    end

    c = calcul_dist(prp)
    
    n = prp["Clients"].nb_clients      # Nombre de clients (usine inclue)
    l = prp["l"]                       # Nombre de périodes - t appartient(0, ..., l)
    m = prp["k"]
    
    @variable(LP, x_dec[1:n, 1:n], Bin)
    @variable(LP, w_dec[1:n] >= 0)
    
    # Objective Function
    
    @objective(LP, Min, sum(sum(c[i, j] * x_dec[i, j] for i in 1:n if i ≠ j) for j in 1:n))
    
    for i in 1:n
        delete(LP, x_dec[i, i])
    end
    
    # Constraint
    @constraint(LP, sum(x_dec[1, j] for j in 2:n) <= m) # (6)
    @constraint(LP, sum(x_dec[i, 1] for i in 2:n) <= m) # (7)
    for i in 2:n
        @constraint(LP, sum(x_dec[i, j] for j in 1:n if j ≠ i) == 1) # (8)
        @constraint(LP, sum(x_dec[j, i] for j in 1:n if j ≠ i) == 1) # (9)
        for j in 2:n
            if j ≠ i
                @constraint(LP, w_dec[i] - w_dec[j] >= prp["Clients"].d[i][t] - (prp["Q"] + prp["Clients"].d[i][t]) * (1 - x_dec[i, j])) # (9)
            end
        end
        @constraint(LP, w_dec[i] <= prp["Q"]) # (9)
    end

    optimize!(LP)
    return value.(x), value.(w)
end