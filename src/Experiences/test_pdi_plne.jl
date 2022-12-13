using Dates
using CSV
using DataFrames
Pkg.add("ProgressBars")

include("../readPRP.jl")
include("../PLNE/PDI_BC.jl")
include("../Heuristic/PDI.jl")


function test_heuristics_csv(nb_files, dir, fixed_files)
    res = ""
    log = string("log/log_", Dates.now(), ".csv")
    files = readdir(dir)[1:nb_files]

    
    df = DataFrame(file = [], 
                PLNE_score = [],
                PLNE_time = [],
                )
    
    i = 0

    for t in ProgressBar(1:nb_files+length(fixed_files))
        if t > nb_files
            file = fixed_files[t-nb_files]
            f = string(dir, "/", file)
        else
            file = files[t]
            f = string(dir, "/", file)
        end
        s = zeros(3)
        time = zeros(3)

        if i == 0
            (c_bp, lbp), t_bp = @timed PDI_H(f, "BP") # JUST TO REMOVE THE FIRST FAKE VALUE
        end

        i += 1
        
        prp = Read_PRP_instance(f)

        (p_dec, y_dec, I_dec, q_dec, x_dec, l, n, cost), t_plne = @timed PDI_BC(prp)
            
        s[1] += cost
        time[1] += t_plne
        
        df_current = DataFrame(file = [file], PLNE_score = [s[1]], PLNE_time = [time[1]])
        append!(df, df_current)
    end
    CSV.write(log, df)
    
    return df
end