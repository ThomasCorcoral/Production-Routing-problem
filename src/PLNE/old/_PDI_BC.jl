using JuMP
using CPLEX
using HiGHS
using Combinatorics

include("PDI1NoBranch.jl")

function PDI_BC(prp, optim="CPLEX", time_limit=3600, gfsec=true, fcc=false, user_cut=false, reverse=false)
    
    LP = PDI1(prp, true, optim)
    
    function lazySep(c_data)
        
        n = prp["n"]+1
        l = prp["l"]
        Q = prp["Q"]
        
        for t in 1:l
            fullH = Set{Int64}()
            sts = Set{Int64}[]
            for i in 1:n
                if !(i in fullH) && (round(callback_value(c_data, variable_by_name(LP, "z_dec[$i,$t]"))) == 1)
                    foundZero = false 
                    current = i
                    h = Set{Int64}(current)
                    #push!(h, current)
                    while true
                        for j in 1:n
                            if current != j
                                """print(j)
                                print(" ")
                                print(i)
                                print(" ")
                                print(current)
                                println()"""
                                if round(callback_value(c_data, variable_by_name(LP, "x_dec[$j,$current,$t]"))) == 1
                                    current = j
                                    break
                                end
                            end
                        end
                        if current == 0 # Pas de cycle
                            foundZero = true
                            break
                        elseif current in h # Cycle déjà trouvé
                            break
                        end
                        push!(h, current) # Nouveau cycle
                    end
                    if !foundZero # ! Sous tour
                        push!(sts, h)
                    end
                    fullH = union(fullH, h) # Ajout du cycle
                end
            end
            
            fullSet = Set(1:n)
            for st in sts
                not_st = setdiff(fullSet, st)
                if fcc
                    #FCCs
                    if length(st)>=1
                        con = @build_constraint(sum(variable_by_name(LP, "x_dec[$i,$j,$t]") for i in not_st, j in st) >= sum(variable_by_name(LP, "q_dec[$i,$t]") for i in st) / Q)
                        MOI.submit(LP, MOI.LazyConstraint(c_data), con)
                    end
                    if reverse
                        #FCCs Reversed Q
                        con = @build_constraint(Q * sum(variable_by_name(LP, "x_dec[$i,$j,$t]") for i in not_st, j in st) >= sum(variable_by_name(LP, "q_dec[$i,$t]") for i in st))
                        MOI.submit(LP, MOI.LazyConstraint(c_data), con)
                    end
                end
                
                if gfsec
                    if length(st)>=2
                        con = @build_constraint(sum( i != j ? variable_by_name(LP, "x_dec[$i,$j,$t]") : 0 for i in st, j in st) <= length(st) - sum( variable_by_name(LP, "q_dec[$i,$t]") for i in st) / Q)
                        MOI.submit(LP, MOI.LazyConstraint(c_data), con)
                    end
                    if reverse
                        if length(st)>=2
                            con = @build_constraint(Q * sum( i != j ? variable_by_name(LP, "x_dec[$i,$j,$t]") : 0 for i in st, j in st) <= sum(Q * variable_by_name(LP, "z_dec[$i,$t]") - variable_by_name(LP, "q_dec[$i,$t]") for i in st))
                            MOI.submit(LP, MOI.LazyConstraint(c_data), con)
                        end
                    end
                end
                
                if user_cut
                    MOI.submit(LP, MOI.UserCut(c_data), con) 
                end
            end
        end
    end

    MOI.set(LP, MOI.LazyConstraintCallback(), lazySep)
    set_time_limit_sec(LP, time_limit) # 1H = 3600
    optimize!(LP)
    return LP
end