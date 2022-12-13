using JuMP
using CPLEX

# Cf. Project subject
function LSP(data, SC)

    n = data["n"]
    l = data["l"]
    M = data["C"]

    LP = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_EPINT" => 1e-15 ))

    @variable(LP, 0 <= p_dec[1:l])
    @variable(LP, y_dec[1:l], Bin)
    @variable(LP, 0 <= I_dec[0:n, 0:l])
    @variable(LP, 0 <= q_dec[1:n, 1:l])
    if SC != 0
        @variable(LP, z_dec[1:n, 1:l], Bin) 
    end
    if SC == 0
        @objective(LP, Min, sum(data["u"] * p_dec[t] + data["f"] * y_dec[t] + sum(data["Clients"].h[i+1] * I_dec[i, t] for i = 1:n) for t in 1:l))
    else
        @objective(LP, Min, sum(data["u"] * p_dec[t] + data["f"] * y_dec[t] + sum(data["Clients"].h[i+1] * I_dec[i, t] + SC[i, t] * z_dec[i, t] for i in 1:n) for t in 1:l))
        for i in 1:n
            for t in 1:l
                @constraint(LP, q_dec[i, t] <= z_dec[i, t] * M)
            end
        end
    end
   for i in 0:n
        @constraint(LP, I_dec[i, 0] == data["Clients"].L0[i+1]) # Init inventory at t = 0 / T
   end

   for t in 1:l
        @constraint(LP, I_dec[0, t-1] + p_dec[t] == sum(q_dec[i, t] for i in 1:n) + I_dec[0, t]) # (1)
        for i in 1:n
            @constraint(LP, I_dec[i, t-1] + q_dec[i, t] == data["Clients"].d[i, t] + I_dec[i, t] ) # (2)
            @constraint(LP, I_dec[i, t-1] + q_dec[i, t] <= data["Clients"].L[i+1] ) # (5)
        end
        @constraint(LP, p_dec[t] <= M * y_dec[t]) # (3)
        @constraint(LP, I_dec[0, t-1] <= data["Clients"].L[1]) # (4) 
        @constraint(LP, sum(q_dec[i, t] for i in 1:n) <= data["Q"] * data["k"]) 
    end
    set_silent(LP)
    optimize!(LP)
    return value.(p_dec), value.(y_dec), value.(I_dec), value.(q_dec)
end