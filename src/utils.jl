using Graphs


function dist(C, i, j, t, mc)
    if(t == 1)
        return floor.(((C.X[i] - C.X[j])^2 + (C.Y[i] - C.Y[j])^2)^(0.5) + 0.5)
    else
        return ((C.X[i] - C.X[j])^2 + (C.Y[i] - C.Y[j])^2)^(0.5) * mc
    end
end


function calcul_dist(p)
    if(p["Type"] == 1)
        mc = 0
    else
        mc = p["mc"]
    end
    c = Array{Float64}(undef, (p["n"]+1, p["n"]+1))
    for i in 1:p["n"]+1
        for j in 1:p["n"]+1            
            c[i, j] = dist(p["Clients"], i, j, p["Type"], mc)
        end
    end
    return c
end

function find_tournee(n,xkt,i,circuit)
	"""
	circuit: un chiffre, le point de depart du circuit
	n : nb de revendeur
	i : le point de depart suivant
	xkt : matrice x de la voiture k, a l instant t
	donne le circuit de la voirture k a l instant t
	"""
	depart = i
	tournee = [circuit]
	findDepot = false
	dejaPasse = Set() #stock les noeuds deja passes
	push!(dejaPasse, circuit)
	#println(xkt)
	while(!findDepot && !(depart in dejaPasse))
		push!(tournee, depart)
		push!(dejaPasse, depart)
		for j = 1:n+1
			if(xkt[depart, j] > 0)
				depart = j
				break
			end
		end
		if(depart == 1)
			if(!(1 in tournee))
				push!(tournee, depart)
			end
			findDepot = true
		end
	end
	return findDepot, tournee
end

function Mt(data,t)
	s=0
	for T in [t : data["l"];]
		for i in [1 : data["n"];]
			s=s+data["d"][i,T]
		end
	end
	return min(data["C"],s)
end

function Mit_til(data,i,t)
	s=0
	for j in [t : data["l"];]
		s=s+data["d"][i,j]
	end
	return min(data["L"][i],data["Q"],s)
end

function detect_sous_tour(n, xkt)
	"""
	n: nb de revendeurs
	xkt : matrice x[:,:,k,t], contient peut etre plusieurs tournees
	detecter un sous tour d une voiture k, a l instant t
	"""
	sousTours=Set()
	dejaVisite=Set()
	trouveDepartDepot=true # =false si trouver au moins un sous tour
	for i = 2:n+1
		if(!(i in dejaVisite) ) #si deja visite, appartient deja a un sous tours
			for j =2:n+1
				if(xkt[i,j]>=1)
					passeDepot,tournee=find_tournee(n,xkt,j,i)
					if(!passeDepot && !(1 in tournee)) #detection d un sous tours
						trouveDepartDepot=false
						push!(sousTours, Set(tournee))
						union!(dejaVisite,tournee)

					end
				end
			end
		end
	end
	return trouveDepartDepot,sousTours
end

function cost_circuit(cout, circuit)
    c = circuit[1]
    cout_total = 0
    for (from, to) in zip(c[1:end-1], c[2:end])
        cout_total += cout[from+1, to+1]
    end
    cout_total += cout[c[end]+1, c[1]+1]
    return cout_total
end

function generate_coords(prp, dict, clients_t)
	coord = [(prp["Clients"].X[i], prp["Clients"].Y[i]) for i in 1:prp["n"]+1]
	coordx_t = [0 for i in 1:length(clients_t)]
	coordy_t = [0 for i in 1:length(clients_t)]
	for i in clients_t
		coordx_t[dict[i]] = coord[i+1][1]
		coordy_t[dict[i]] = coord[i+1][2]
	end
	return coordx_t, coordy_t
end