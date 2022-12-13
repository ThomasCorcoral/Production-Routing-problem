using Dates
using CSV
using DataFrames
using ProgressBars

# include("../utils.jl")
include("../Heuristic/PDI.jl")


function test_heuristics_csv(nb_files, dir, fixed_files)
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

        (c_bp, lbp), t_bp = @timed PDI_H(f, "BP")
        (c_cw, lcw), t_cw = @timed PDI_H(f, "CW")
        (c_sect, lsect), t_sect = @timed PDI_H(f, "SECT")
            
        s[1] += c_bp
        s[2] += c_cw
        s[3] += c_sect
        time[1] += t_bp
        time[2] += t_cw
        time[3] += t_sect
        
        df_current = DataFrame(file = [file], BP_score = [s[1]], CW_score = [s[2]], SECT_score = [s[3]], BP_time = [time[1]], CW_time = [time[2]], SECT_time = [time[3]])
        append!(df, df_current)
    end
    CSV.write(log, df)
    
    return df
end