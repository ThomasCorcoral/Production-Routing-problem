using JuMP
using CPLEX
using CPUTime
using Graphs
import HiGHS

include("../utils.jl")


function PDI1(prp, gfsec, optim="CPLEX")
    if optim == "HiGHS"
        LP = Model(HiGHS.Optimizer)
    elseif optim == "CPLEX"
        LP = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_EPINT" => 1e-15 ))
    else
        return -1
    end

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
    
    # Constraint
    
    for t in 1:l
        if t>1
            @constraint(LP, I_dec[1, t-1] + p_dec[t] == sum(q_dec[i, t] for i in 2:n) + I_dec[1, t]) # (2)
        else
            @constraint(LP, prp["Clients"].L0[1] + p_dec[t] == sum(q_dec[i, t] for i in 2:n) + I_dec[1, t]) # (2)
        end
        for i in 1:n
            @constraint(LP, x_dec[i, i, t] == 0)
            if i > 1
                if t>1
                    @constraint(LP, I_dec[i, t-1] + q_dec[i, t] == prp["Clients"].d[i][t] + I_dec[i, t]) # (3)
                    @constraint(LP, I_dec[i, t-1] + q_dec[i, t] <= prp["Clients"].L[i]) # (6)
                else
                    @constraint(LP, prp["Clients"].L0[i] + q_dec[i, t] == prp["Clients"].d[i][t] + I_dec[i, t]) # (3)
                    @constraint(LP, prp["Clients"].L0[i] + q_dec[i, t] <= prp["Clients"].L[i]) # (6)
                end
                Mit = min(prp["Clients"].L[i], prp["Q"], sum(prp["Clients"].d[i][j] for j in t:l))
                @constraint(LP, q_dec[i, t] <= Mit * z_dec[i, t]) # (7)
                @constraint(LP, sum(x_dec[i, j, t] for j in 1:n) == z_dec[i, t]) # 8
                if(!gfsec)
                    @constraint(LP, 0 <= w_dec[i, t]) # (12)
                    @constraint(LP, w_dec[i, t] <= prp["Q"] * z_dec[i, t]) # (12)
                end
                set_binary(z_dec[i,t]) # (15)
            end
            @constraint(LP, sum(x_dec[i, j, t] for j in 1:n) + sum(x_dec[j, i, t] for j in 1:n) == 2*z_dec[i, t]) # (9)
            # (11) START
            Mit = min(prp["Clients"].L[i], prp["Q"], sum(prp["Clients"].d[i][j] for j in t:l))
            if(!gfsec)
                for j in 1:n
                   if i ≠ j
                       @constraint(LP, w_dec[i, t] - w_dec[j, t] >= q_dec[i, t] - Mit * (1 - x_dec[i, j, t]))
                   end
                end
            end
            # (11) END
            @constraint(LP, 0 <= I_dec[i, t]) # (13)
            @constraint(LP, 0 <= q_dec[i, t]) # (13)
            # (14) Defined as a binary variable
        end
        Mt = min(prp["C"], sum(sum(prp["Clients"].d[i][j] for i in 2:n) for j in t:l))
        @constraint(LP, p_dec[t] <= Mt * y_dec[t]) # (4)
        @constraint(LP, I_dec[1, t] <= prp["Clients"].L[1]) # (5)
        @constraint(LP, z_dec[1, t] <= prp["k"]) # (10)
        @constraint(LP, 0 <= p_dec[t]) # (13)
        @constraint(LP, 0 <= z_dec[1, t]) # (16)
        set_integer(z_dec[1, t]) # (16)
    end
    return LP
end
