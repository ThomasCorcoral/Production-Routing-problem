using JuMP
using CPLEX
#using CPUTime
#using Graphs
#import HiGHS

include("../utils.jl")

function PDI1(prp)
    #LP = Model(HiGHS.Optimizer)
    LP = Model(CPLEX.Optimizer)
    
    c = calcul_dist(prp)
    
    n = prp["Clients"].nb_clients      # Nombre de clients (usine inclue)
    l = prp["l"]                       # Nombre de périodes - t appartient(0, ..., l)
    
    @variable(LP, y_dec[1:l], Bin)
    @variable(LP, p_dec[1:l])
    @variable(LP, I_dec[1:n, 1:l], Int)
    @variable(LP, x_dec[1:n, 1:n, 1:l], Bin)
    @variable(LP, q_dec[1:n, 1:l])
    @variable(LP, z_dec[1:n, 1:l])
    @variable(LP, w_dec[1:n, 1:l])
    
    # Objective Function (1)
    
    @objective(LP, Min, sum(prp["u"] * p_dec[t] + 
            prp["f"] * y_dec[t] + 
            sum(prp["Clients"].h[i] * I_dec[i, t] for i in 1:n) + 
            sum(sum(c[i, j] * x_dec[i, j, t] for i in 1:n if i ≠ j) for j in 1:n) for t in 1:l))
    
    # Constraint (2)
    
    for t in 1:l
        if t>1
            @constraint(LP, I_dec[1, t-1] + p_dec[t] == sum(q_dec[i, t] for i in 2:n) + I_dec[1, t])
        else
            @constraint(LP, prp["Clients"].L0[1] + p_dec[t] == sum(q_dec[i, t] for i in 2:n) + I_dec[1, t])
        end
    end
    
    # Constraint (3)
    
    for t in 1:l
        for i in 2:n
            if t>1
                @constraint(LP, I_dec[i, t-1] + q_dec[i, t] == prp["Clients"].d[i][t] + I_dec[i, t])
            else
                @constraint(LP, prp["Clients"].L0[i] + q_dec[i, t] == prp["Clients"].d[i][t] + I_dec[i, t])
            end
        end
    end
    
    # Constraint (4)
    
    for t in 1:l
        Mt = min(prp["C"], sum(sum(prp["Clients"].d[i][j] for i in 2:n) for j in t:l))
        @constraint(LP, p_dec[t] <= Mt * y_dec[t]) 
    end
    
    # Constraint (5)
    
    for t in 1:l
        @constraint(LP, I_dec[1, t] <= prp["Clients"].L[1])
    end
    
    # Constraint (6)
    
    for t in 1:l
        for i in 2:n
            if t > 1
                @constraint(LP, I_dec[i, t-1] + q_dec[i, t] <= prp["Clients"].L[i])
            else
                @constraint(LP, prp["Clients"].L0[i] + q_dec[i, t] <= prp["Clients"].L[i])
            end
        end
    end
    
    # Constraint (7)
    
    for t in 1:l
        for i in 2:n
            Mit = min(prp["Clients"].L[i], prp["Q"], sum(prp["Clients"].d[i][j] for j in t:l))
            @constraint(LP, q_dec[i, t] <= Mit * z_dec[i, t])
        end
    end
    
    # Constraint (8)
    
    for t in 1:l
        for i in 2:n
            @constraint(LP, sum(x_dec[i, j, t] for j in 1:n) == z_dec[i, t])
        end
    end
    
    # Constraint (9)
    
    for t in 1:l
        for i in 1:n
            @constraint(LP, sum(x_dec[i, j, t] for j in 1:n) + sum(x_dec[j, i, t] for j in 1:n) == 2*z_dec[i, t])
        end
    end
    
    # Constraint (10)
    
    for t in 1:l
        @constraint(LP, z_dec[1, t] <= prp["k"])
    end
    
    # Constraint (11)
    
    for t in 1:l
        for i in 1:n
            Mit = min(prp["Clients"].L[i], prp["Q"], sum(prp["Clients"].d[i][j] for j in t:l))
            for j in 1:n
                if i ≠ j
                    @constraint(LP, w_dec[i, t] - w_dec[j, t] >= q_dec[i, t] - Mit * (1 - x_dec[i, j, t]))
                end
            end
        end
    end
    
    # Constraint (12)
    
    for t in 1:l
        for i in 2:n
            @constraint(LP, 0 <= w_dec[i, t])
            @constraint(LP, w_dec[i, t] <= prp["Q"] * z_dec[i, t])
        end
    end
    
    # Constraint (13)
    
    for t in 1:l
        @constraint(LP, 0 <= p_dec[t])
        for i in 1:n
            @constraint(LP, 0 <= I_dec[i, t])
            @constraint(LP, 0 <= q_dec[i, t])
        end
    end
    
    # Constraint (14)
    
    # Defined as a binary variable
    
    # Constraint (15)
    
    for t in 1:l
        for i in 2:n
            #@constraint(LP, z_dec[i, t]
            set_binary(z_dec[i,t])
        end
    end
    
    # Constraint (16)
    
    for t in 1:l
        @constraint(LP, 0 <= z_dec[1, t]) 
        set_integer(z_dec[1, t])
    end
    
    # Remove x[i,i,t]
    for t in 1:l
        for i in 1:n
            @constraint(LP, x_dec[i, i, t] == 0)
        end
    end
    # End
    
    # print(LP)
    # println()
    # optimize!(LP)
   
    # println(solution_summary(LP, verbose=true))
    
    return LP
end