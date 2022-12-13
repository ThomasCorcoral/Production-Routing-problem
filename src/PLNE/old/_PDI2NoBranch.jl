using JuMP
using CPLEX
#using CPUTime
#using Graphs
#import HiGHS

include("../utils.jl")

function PDI2(prp)
    #LP = Model(HiGHS.Optimizer)
    LP = Model(CPLEX.Optimizer)

    """n = params["n"] #nombre de revendeur (exclu le fournisseur donc)
    m = params["k"] #nb de véhicules
    l=params["l"] #nombre de pas de temps
    Q = params["Q"] #capacité d'un vahicule
    C=params["C"] #production max du fournisseur à un certain pas de temps
    u=params["u"] #coût de production d'une unité
    f=params["f"] #coût de setup
    a=1"""
    
    c = calcul_dist(prp)
    #n = prp["n"]                        # Nombre de clients (usine inclue)
    n = prp["Clients"].nb_clients    # Nombre de clients (usine inclue)
    l = prp["l"]                        # Nombre de périodes - t appartient(0, ..., l)
    m = prp["k"]
    m = 200

    # Variables

    """@variable(LP, p_dec[t=1:l] >= 0)
    @variable(LP, I_dec[0:n, 0:l] >= 0)
    @variable(LP, q_dec[1:n, k=1:m, 1:l] >= 0)
    @variable(LP, y_dec[1:l], Bin)
    @variable(LP, x_dec[i=0:n, j=0:n, k=1:m, t=1:l], Bin)
    @variable(LP, z_dec[i=0:n, k=1:m, t=1:l], Bin)"""

    @variable(LP, y_dec[1:l], Bin)
    @variable(LP, p_dec[1:l] >= 0)
    @variable(LP, I_dec[1:n, 1:l] >= 0, Int)
    @variable(LP, x_dec[1:n, 1:n, 1:m, 1:l], Bin)
    @variable(LP, q_dec[1:n, 1:m, 1:l] >= 0)
    @variable(LP, z_dec[1:n, 1:m, 1:l], Bin)

    # Objectif

    @objective(LP, Min, sum(prp["u"] * p_dec[t] + 
            prp["f"] * y_dec[t] + 
            sum(prp["Clients"].h[i] * I_dec[i, t] for i in 1:n) + 
            sum(sum(c[i, j] * sum(x_dec[i, j, k, t] for k in 1:m) for i in 1:n if i ≠ j) for j in 1:n) for t in 1:l))

    # Contrainte (21)
    
    for t in 1:l
        if t>1
            @constraint(LP, I_dec[1, t-1] + p_dec[t] == sum(sum(q_dec[i, k, t] for k in 1:m) for i in 2:n) + I_dec[1, t])
        else
            @constraint(LP, prp["Clients"].L0[1] + p_dec[t] == sum(sum(q_dec[i, k, t] for k in 1:m) for i in 2:n) + I_dec[1, t])
        end
    end
    
    # Contrainte (22)
    
    for t in 1:l
        for i in 2:n
            if t>1
                @constraint(LP, I_dec[i, t-1] + sum(q_dec[i, k, t] for k in 1:m) == prp["Clients"].d[i][t] + I_dec[i,t])
            else
                @constraint(LP, prp["Clients"].L0[1] + sum(q_dec[i, k, t] for k in 1:m) == prp["Clients"].d[i][t] + I_dec[i,t])
            end
        end
    end
    
    # Contrainte (23)
    
    for t in 1:l
        Mt = min(prp["C"], sum(sum(prp["Clients"].d[i][j] for i in 2:n) for j in t:l))
        @constraint(LP, p_dec[t] <= Mt * y_dec[t]) 
    end
    
    # Contrainte (24)
    
    for t in 1:l
        @constraint(LP, I_dec[1, t] <= prp["Clients"].L[1])
    end
    
    # Contrainte (25)
    
    for t in 1:l
        for i in 2:n
            if t > 1
                @constraint(LP, I_dec[i, t-1] + sum(q_dec[i, k, t] for k in 1:m) <= prp["Clients"].L[i])
            else
                @constraint(LP, prp["Clients"].L0[i] + sum(q_dec[i, k, t] for k in 1:m) <= prp["Clients"].L[i])
            end
        end
    end
    
    # Contrainte (26)
    
    for t in 1:l
        for i in 2:n
            Mit = min(prp["Clients"].L[i], prp["Q"], sum(prp["Clients"].d[i][j] for j in t:l))
            for k in 1:m
                @constraint(LP, q_dec[i, k, t] <= Mit * z_dec[i, k, t])
            end
        end
    end
    
    # Contrainte (27)
    
    for t in 1:l
        for i in 2:n
            @constraint(LP, sum(z_dec[i, k, t] for k in 1:m) <= 1)
        end
    end
    
    # Contrainte (28)
    
    for t in 1:l
        for i in 1:n
            for k in 1:m
                @constraint(LP, sum(x_dec[i, j, k, t] for j in 1:n) + sum(x_dec[j, i, k, t] for j in 1:n) == 2*z_dec[i, k, t])
            end
        end
    end
    
    # Contrainte (29)
    for t in 1:l
        for i in 2:n
            for k in 1:m
                for S in powerset(2:n, 2, n)
                    @constraint(LP, sum(sum(x_dec[i, j, k, t] for j in S) for i in S) <= length(S)-1) #(29)
                end
            end
        end
    end
    
    # Contrainte (30)
    
    for t in 1:l
        for k in 1:m
            @constraint(LP, sum(q[i, k, t] for i in 1:n) <= prp["Q"] * z[0, k, t]) #(30)
        end
    end
    # End
    
    # print(LP)
    # println()
    # optimize!(LP)
   
    # println(solution_summary(LP, verbose=true))
    
    return LP
end