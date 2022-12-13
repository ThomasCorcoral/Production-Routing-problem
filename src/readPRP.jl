


struct Clients
    nb_clients
    X
    Y
    h
    L
    L0
    d
end


function Read_PRP_instance(filename, plne=true)
    hashmap = Dict()
    open(filename) do f
        
        demands = false
        listing_active = false
        current_value = 0
        
        nbc = 0
        X = Array{Float64}(undef, 0)
        Y = Array{Float64}(undef, 0)
        h = Array{Float64}(undef, 0)
        L = Array{Float64}(undef, 0)
        L0 = Array{Float64}(undef, 0)
        d = Array{Array{Float64}}(undef, 0)
        
        for (i, line) in enumerate(eachline(f))
            x = split(line, " ")
            deleteat!(x, findall(e->e=="", x))
            current_item = x[1]
            
            if(current_item == "0" || listing_active) # Liste de clients / demandes
                listing_active = true
                if(current_item == "d")
                    push!(d, zeros(hashmap["l"]+1))
                    demands = true
                elseif(demands) # Enregistrement des demandes
                    push!(d, append!([0],parse.(Float64, x)[2:end]))
                else # enregistrement des clients
                    nbc += 1
                    push!(X, parse(Float64, x[2]))
                    push!(Y, parse(Float64, x[3]))
                    push!(h, parse(Float64, x[6]))
                    push!(L, parse(Float64, x[8]))
                    push!(L0, parse(Float64, x[10]))
                end
            else # variables générales du problème
                try
                    current_value = parse(Int64, x[2])
                catch err
                    current_value = 10000000000
                end
                hashmap[current_item] = current_value
            end
        end
        if !plne
            d = reduce(hcat,d)'[2:end, 2:end]
        end
        C = Clients(nbc, X, Y, h, L, L0, d)
        hashmap["Clients"] = C
    end
    return hashmap
end