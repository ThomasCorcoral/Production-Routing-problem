using JuMP
using CPLEX
using DelimitedFiles
using Dates

include("VRP.jl")

include("../utils.jl")
include("../draw.jl")
include("../PLNE/LSP.jl")
include("../readPRP.jl")

# PDI function with heuristic
function PDI_H(file, h)
    """
    # Arguments
    - 'file::String': file path with the prp instance.
    - 'h::String': The heuristic {"BP", "CW", "SECT"}.
    """

    prp = Read_PRP_instance(file, false)

    itMax = 3
    n = prp["n"]
    l = prp["l"]
    
    c = calcul_dist(prp)

    SC = Array{Float64}(undef,n,l)
    for i=1:n
        for t=1:l
            SC[i,t] = c[1,i+1] + c[i+1,1]
        end
    end

    bestSol = []
    best_cost = 0
    it = 0
    while it < itMax
        # LSP #
        p_dec, y_dec, I_dec, q_dec = LSP(prp, SC)
        CoutSolcourante = sum(prp["u"] * p_dec[t] + prp["f"] * y_dec[t] + sum(prp["Clients"].h[i] * I_dec[i-1, t] for i = 1:n+1) for t = 1:l )

        # VRP #
        tours, CoutSolcourante = process_VRP(prp, q_dec, h, CoutSolcourante)

        if it == 0 || best_cost > CoutSolcourante
            bestSol = [p_dec, y_dec, I_dec, q_dec, tours]
            best_cost = CoutSolcourante
        end
        for t_tournee in 1:length(tours)
            elems_t = [] # Delivery list
            for cycle in tours[t_tournee] 
                for ind_elem in 1:length(cycle) 
                    push!(elems_t, cycle[ind_elem])
                    if ind_elem == 1 || cycle[ind_elem] == 0
                        continue
                    elseif ind_elem == length(cycle)
                        SC[cycle[ind_elem], t_tournee] = c[cycle[ind_elem-1]+1, cycle[ind_elem]+1] + c[cycle[ind_elem]+1, cycle[1]+1] - c[cycle[ind_elem-1]+1, cycle[1]+1]
                    else
                        SC[cycle[ind_elem], t_tournee] = c[cycle[ind_elem-1]+1, cycle[ind_elem]+1] + c[cycle[ind_elem]+1, cycle[ind_elem+1]+1] - c[cycle[ind_elem-1]+1, cycle[ind_elem+1]+1]
                    end 
                end
            end        
            for i in 1:n
                if !(i in elems_t)
                    SC = update_SC(tours, t_tournee, SC, i, c)
                end
            end
        end
        it += 1
    end

    p, y, I, q, tours = bestSol
    date = string(Dates.now())
    f = collect(eachsplit(collect(eachsplit(file, "/"))[end], "."))[1]
    for t in 1:length(tours)
        if tours[t] != [] # If there is a tour at this period
            pdi_to_png(prp, q, t, tours[t], f, date)
        end
    end
    return best_cost, l
end

# function to process the Vehicule Routing Problem
# return the tours and the cost
# Used in function PDI_H
function process_VRP(prp, q_dec, h, cost)
    tours = []
    final_cost = cost
    for t=1:prp["l"]
        # START VRP
        if sum(q_dec[i,t] for i in 1:prp["n"]) > 0
            cycles_t = VRP_iteratif(prp, q_dec, t, h)
            push!(tours, cycles_t)
            final_cost += get_vrp_cost(prp, cycles_t)
        else # NO TOUR
            push!(tours,[]) 
        end
    end
    return tours, final_cost
end

# Trivial function - update the SC and return it
# Used in function PDI_H
function update_SC(tours, t, SC, i, c)
    if tours[t] == []              
        SC[i, t]= c[1,i+1] + c[i+1,1]
    else
        cycle = tours[t][1]
        if length(cycle) == 1
            _min = c[cycle[1]+1, i+1] + c[i+1, cycle[1]+1]
        else
            _min = c[cycle[1]+1, i+1] + c[i+1, cycle[2]+1] - c[cycle[1]+1, cycle[2]+1]
        end
        for cycle in tours[t] 
            for ii in 1:length(cycle)
                if ii == length(cycle)
                    _cur = c[cycle[ii]+1, i+1] + c[i+1, cycle[1]+1] - c[cycle[ii]+1, cycle[1]+1]
                elseif ii != length(cycle)
                    _cur = c[cycle[ii]+1, i+1] + c[i+1, cycle[ii+1]+1] - c[cycle[ii]+1, cycle[ii+1]+1]
                end
                if _cur < _min
                    _min = _cur 
                end
            end
        end
        SC[i, t]= _min
    end
    return SC
end

function get_vrp_cost(data, cycles)
    cost = calcul_dist(data)
    res = 0
    for c in cycles
        res += get_cycle_cost(cost, c)	
    end
    return res
end