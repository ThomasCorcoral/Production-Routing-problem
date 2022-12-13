using Dates
using CSV
using DataFrames

include("../utils.jl")
include("../readPRP.jl")
include("../Heuristic/binPacking.jl")
include("../Heuristic/clarkWright.jl")
include("../Heuristic/sectorielle.jl")

function test_heuristics_csv(nb_files, dir)
    plne=false
    res = ""
    log = string("log/log_", Dates.now(), ".csv")
    files = readdir(dir)[1:nb_files]
    
    df = DataFrame(file = [], 
                BP_score = [],
                CW_score = [],
                SECT_score = [],
                BP_time = [],
                CW_time = [],
                SECT_time = [],
                )
    
    for file in files
        f = string(dir, "/", file)
        h = Read_PRP_instance(f)
        s = zeros(4)
        time = zeros(4)
        if plne
            pl, t_pdi = @timed PDI_BC(h)
            x = all_variables(pl)
            x = x[end-(h["n"]+1)*(h["n"]+1)*h["l"]+1:end]
        end
        
        costs = calcul_dist(h)
        for t in 1:h["l"]
            sect, t_sect = @timed Sectorielle(h, t)
            cw, t_cw = @timed clark_wright(h, h["Clients"].d, costs, t)
            bp, t_bp = @timed binPacking(h, t)
            
            if plne
                c_i = 1
                c_j = 1
                c_plne = 0
                for i in (h["n"]+1)*(h["n"]+1)*(t-1)+1 : (h["n"]+1)*(h["n"]+1)*(t)
                    c_plne += value.(x[i]) * costs[c_i, c_j]
                    c_i = (c_i + 1) % (h["n"]+2)
                    if c_i == 0
                        c_i += 1
                        c_j += 1
                    end
                end
                s[4] += c_plne
                time[4] += t_pdi
            end
            c_bp = cost_circuit(costs, bp)
            c_cw = cost_circuit(costs, cw)
            c_sect = cost_circuit(costs, sect)
            
            s[1] += c_bp
            s[2] += c_cw
            s[3] += c_sect
            time[1] += t_bp
            time[2] += t_cw
            time[3] += t_sect
        end
        
        df_current = DataFrame(file = [file], BP_score = [s[1]], CW_score = [s[2]], SECT_score = [s[3]], BP_time = [time[1]], CW_time = [time[2]], SECT_time = [time[3]])
        append!(df, df_current)
    end
    CSV.write(log, df)
    
    return df
end