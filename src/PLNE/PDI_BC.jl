using IterTools
using JuMP
using CPLEX
using Combinatorics

include("PDI1NoBranch.jl")
include("../utils.jl")

function PDI_BC(prp)

    c = calcul_dist(prp)
    LP = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_EPINT" => 1e-15 ))


    n = prp["n"]                       # Nombre de clients (usine exclue)
    l = prp["l"]                       # Nombre de périodes - t appartient(0, ..., l)
    m = prp["k"]


    @variable(LP, 0 <= p_dec[1:l])
    @variable(LP, y_dec[1:l], Bin)
    @variable(LP, z_dec[1:n, 1:l], Bin)
    @variable(LP, 0 <= z0_dec[1:l] <= m, Int)
    @variable(LP, 0 <= w_dec[1:n, 1:l] )

    @variable(LP, 0 <= I_dec[1:n+1, 1:l])
    @variable(LP, 0 <= q_dec[1:n+1, 1:l], Int)
    @variable(LP, x_dec[1:n+1, 1:n+1, 1:l], Bin)

    @objective(LP, Min, sum(prp["u"]*p_dec + 
                        prp["f"]*y_dec + 
                        sum(prp["Clients"].h .* I_dec, dims=1)[1,:] + 
                        sum(sum(c .* x_dec, dims=1), dims=2)[1, 1, :])) 

    #contraintes
    for t in [1:l;]
        for i =1:n+1
            @constraint(LP,x_dec[i, i, t]==0)
        end
        if(t==1)
            @constraint(LP, prp["Clients"].L0[1] + p_dec[t] == sum(q_dec[:,t]) - q_dec[1, t] + I_dec[1, t] ) # (2)
            for i in [2:n+1;] 
                @constraint(LP, prp["Clients"].L0[i] + q_dec[i, t] == prp["Clients"].d[i][t+1] + I_dec[i, t]) # (3)
                @constraint(LP, prp["Clients"].L0[i] + q_dec[i, t] <= prp["Clients"].L[i]) # (6)
            end
        else 
            @constraint(LP, I_dec[1, t-1] + p_dec[t] == sum(q_dec[:, t]) - q_dec[1, t] + I_dec[1, t]) # (2)
        end
        Mt = min(prp["C"], sum(sum(prp["Clients"].d[i][j+1] for i in 2:n+1) for j in t:l))
        @constraint(LP, p_dec[t] <= Mt * y_dec[t]) # (4)
        @constraint(LP, I_dec[1, t] <= prp["Clients"].L[1]) # (5)
        @constraint(LP, z0_dec[t] <= m) # (10) 
        for i in [1:n+1;]
            Mit = min(prp["Clients"].L[i], prp["Q"], sum(prp["Clients"].d[i][j+1] for j in t:l))
            if(i==1)
                @constraint(LP, sum(x_dec[i, k, t] for k=1:n+1 if i!=k)+ sum(x_dec[k, i, t] for k=1:n+1 if k!=i) == 2 * z0_dec[t]) # (9)
            elseif(i>1)
                if (t>1)
                    @constraint(LP, I_dec[i, t-1] + q_dec[i ,t] == prp["Clients"].d[i][t+1] + I_dec[i, t]) # (3)
                    @constraint(LP, I_dec[i, t-1] + q_dec[i, t] <= prp["Clients"].L[i]) # (6)
                end
                @constraint(LP, q_dec[i, t] <= Mit * z_dec[i-1, t]) # (7)
                @constraint(LP, sum(x_dec[i, k, t] for k = 1:n+1 if(k!=i)) == z_dec[i-1, t]) # (8)
                @constraint(LP, sum(x_dec[i, k, t] for k = 1:n+1 if k!=i) + sum(x_dec[k, i, t] for k =1:n+1 if k!=i) == 2 * z_dec[i-1, t]) # (9)
            end
        end
    end

    function lazySep_Violated(cb_data)
        for t = 1:l
            xsep = zeros((n+1, n+1))
            for i = 1:n+1
                for j = 1:n+1
                    if(round(callback_value(cb_data, x_dec[i, j, t])) == 1)
                        xsep[i, j] = 1
                    end
                end
            end
            trouveDepartDepot, sousTours = detect_sous_tour(n, xsep)
            if(!trouveDepartDepot)
                for tournee in sousTours
                    if(length(tournee) >= 2)
                        con = @build_constraint(sum(x_dec[i, j, t] for i in tournee for j in tournee if i != j) * prp["Q"] <= sum(prp["Q"] * z_dec[i-1, t] - q_dec[i, t] for i in tournee))
                        MOI.submit(LP, MOI.LazyConstraint(cb_data), con)
                        sbar=setdiff(Set(1:n+1), tournee)
                        con = @build_constraint(sum(x_dec[i, j, t] for i in sbar for j in tournee) >= sum(q_dec[i, t] for i in tournee) / prp["Q"]) 
                        MOI.submit(LP, MOI.LazyConstraint(cb_data), con)

                    end
                end
            end
        end
    end
    
    set_silent(LP)
    set_time_limit_sec(LP,600) # 10 minutes
    MOI.set(LP, MOI.LazyConstraintCallback(),lazySep_Violated)
    optimize!(LP)
    # println(solution_summary(LP, verbose=true))
    return value.(p_dec), value.(y_dec), value.(I_dec), value.(q_dec), value.(x_dec), l, n, objective_value(LP)
end

