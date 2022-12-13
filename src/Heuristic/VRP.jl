using Base.Iterators
using LinearAlgebra

include("binPacking.jl")
include("clarkWright.jl")
include("sectorielle.jl")

function switch_heuristic(prp, q_dec, t, heuristic)
	if heuristic == "BP"
		return BP(prp, q_dec, t)
	elseif heuristic == "CW"
		return CW(prp, q_dec, t)
	else
		return SECT(prp, q_dec, t)
	end
end

function get_append_cost(cost, lot, i) 
	cycle = copy(lot)
	push!(cycle, i)
	c_min, i_min = get_cycle_cost(cost, cycle), length(lot)+1
	for ii in 2:length(lot)
		cycle = copy(lot)
		insert!(cycle, ii, i)  
		cur = get_cycle_cost(cost, cycle)
		if c_min > cur
			c_min, i_min = cur, ii
		end
	end
	return c_min, i_min 
end

function get_cycle_cost(cost, cycle)
    res = 0
	for (from, to) in zip(cycle[1:end-1], cycle[2:end])
        res += cost[from+1, to+1]
	end
	res += cost[cycle[end]+1, cycle[1]+1] # DO NOT FORGET TO GET BACK TO FACTORY
	return res
end
 
function TSP(lots, co)
	new_lots = []
	for lot in lots
		new_lot = [0]
		filter!(e->e!=0, lot)
		while length(lot) > 1
			c_min, e_min = co[new_lot[end]+1, lot[1]+1], lot[1]
			for e in lot
				if c_min > co[new_lot[end]+1, e+1]
					c_min, e_min = co[new_lot[end]+1, e+1], e
				end
			end
			push!(new_lot, e_min)
			filter!(e->e!=e_min, lot)
		end
		if length(lot) == 1
			push!(new_lot, lot[1])
		end
		push!(new_lots, new_lot)
	end
	return new_lots
end

function split_routing(lots, k)
	if length(lots) < k 
		l = length(lots)
		for i in 1:k-l
			big = 1
			for ii in 2:length(lots)
				if length(lots[ii]) > length(lots[big])
					big = ii
				end
			end
			cur = lots[big]
			tab1 = cur[1:floor(Int, length(cur) / 2)] 
			tab2 = cur[floor(Int, length(cur) / 2) + 1:length(cur)]
			lots[big] = tab1
			insert!(tab2, 1, 0)
			push!(lots, tab2)
		end
	end
	return lots
end

function vehicles_splitting(lots, k, q_dec, t, Q, cost)
	len = length(lots)
	for i in 1:len-k
		no_change, fail = true, 0
		while no_change
			if fail == 0
				i1, i2 = 1, 1
				for ii in 2:length(lots)
					if length(lots[ii]) < length(lots[i1])
						if length(lots[i2]) > length(lots[i1])
							i2 = i1
						end
						i1 = ii
					end

					if length(lots[ii]) > length(lots[i1]) && length(lots[ii]) < length(lots[i2])
						i2 = ii
					end
				end
				_new = 0
				for e in lots[i1]
					if e != 0
						_new += q_dec[e, t]
					end
				end

				for e in lots[i2]
					if e != 0
						_new += q_dec[e, t]
					end
				end
			
				if _new <= Q && i1 != i2
					# TSP
					lot2 = lots[i2]
					new_lot = [lots[i1] ; lot2[2:length(lot2)]]
					new_lot_tsp = TSP([new_lot], cost)[1]
					lots[i1] =  new_lot_tsp
					deleteat!(lots, i2)
					no_change = false
				end
				fail += 1
			else
				for i1 in 1:length(lots)-1
					if no_change == false
						break
					end
					for i2 in i1+1:length(lots)
						if no_change == false
							break
						end
						_new = 0
						for e in lots[i1]
							if e != 0
								_new += q_dec[e, t]
							end
						end
						for e in lots[i2]
							if e != 0
								_new += q_dec[e, t]
							end
						end
						if _new <= Q
							# TSP ON NEW BOX
							lot2 = lots[i2]
							new_lot = [lots[i1] ; lot2[2:length(lot2)]]
							new_lot_tsp = TSP([new_lot], cost)[1]
							lots[i1] =  new_lot_tsp
							deleteat!(lots, i2) # RM OLD
							no_change = false
						end
						# UNFEASIBLE
						if i1 == length(lots)-1 && i2 == length(lots) && no_change
							return -1
						end
					end
				end
			end
		end
	end
	return lots
end

function lots_ordering(lots, cc, cost, q_dec, t, Q)
	maxit, n_it, change = 5, 0, true

	while n_it < maxit && change
		change = false
		for client in cc
			lot_client = []
			indice_lot_client = -1
			for indice_lot in 1:length(lots) 
				if client in lots[indice_lot]
					lot_client = lots[indice_lot]
					indice_lot_client = indice_lot
					break
				end
			end
			cost_lot_client = get_cycle_cost(cost, lot_client)

			lot_client_sans = copy(lot_client)
			filter!(e->e!=client, lot_client_sans)
			cost_lot_client_sans = get_cycle_cost(cost, lot_client_sans)

			cost_meilleur_change = 0
			indice_meilleure_lot = indice_lot_client
			i_min_best = 0

			for indice_lot in 1:length(lots) 
				if indice_lot != indice_lot_client

					min_cost, min_indice = get_append_cost(cost, lots[indice_lot], client)
					ancien_cost = cost_lot_client + get_cycle_cost(cost, lots[indice_lot])
					nouveau_cost = cost_lot_client_sans + min_cost
					if ancien_cost - nouveau_cost > cost_meilleur_change
						nouvelle_q_dec = 0
						for elemq in lots[indice_lot]
							if elemq != 0
								nouvelle_q_dec += q_dec[elemq, t]
							end
						end
						nouvelle_q_dec += q_dec[client, t]
					
						if nouvelle_q_dec <= Q
							cost_meilleur_change = ancien_cost - nouveau_cost
							indice_meilleure_lot = indice_lot
							i_min_best = min_indice
							change = true 
						end
					end
				end
			end
			if indice_meilleure_lot != indice_lot_client
				new_best = copy(lots[indice_meilleure_lot])
				insert!(new_best, i_min_best, client)

				nouv_lot_sans_client = copy(lot_client)
				filter!(e->e!=client, nouv_lot_sans_client)
				nouv_lots = []
				for indice_lot in 1:length(lots)
					if indice_lot == indice_meilleure_lot  || indice_lot == indice_lot_client
						continue
					end
					push!(nouv_lots, lots[indice_lot])
				end
				lots = nouv_lots
				push!(lots, new_best)
				push!(lots, nouv_lot_sans_client)
			end
		end
        n_it += 1
	end
	return lots
end

function TSP2(lots, cost)
	max_TSP2 = 5 # MAX - n^2
	n_it = 0
	change = [true for i in 1:length(lots)] 
    while n_it < max_TSP2 && maximum(change) # STOP WHEN NO MORE CHANGE OR WHEN MAX ITERATION NUMBER REACH
		for ii in 1:length(lots)
			if change[ii] == false 
				continue
			end
			change[ii] = false
			if length(lots[ii]) > 3
				for i1 in 1:length(lots[ii])-2
					for i2 in i1+2:length(lots[ii])
						if change[ii] # CHECK IF CHANGE ALREADY DONE IN THIS ITERATION
							break
						end
						if i1 != i2 && i1+1 != i2 && i1-1 != i2
							if i1 == 1 && i2 == length(lots[ii])
								continue
							end
							if i2 == length(lots[ii])
								i21 = 1
							else
								i21 = i2 + 1
							end
							# 2-OPT
							if cost[lots[ii][i1]+1, lots[ii][i1+1]+1] + cost[lots[ii][i2]+1, lots[ii][i21]+1] > cost[lots[ii][i1]+1, lots[ii][i2]+1] + cost[lots[ii][i1+1]+1, lots[ii][i21]+1]  							
								neighbour = copy(lots[ii][1:i1])
								# ADD ALL NEIGHBOURS
								for j in i2:-1:i1+1
									push!(neighbour, lots[ii][j])
								end
								# i2 REACH END OF LOTS
								if i2 != length(lots[ii])
									for j in i21:length(lots[ii])
										push!(neighbour, lots[ii][j])
									end 
								end
								lots[ii], change[ii] = neighbour, true
							end
						end
					end
				end
			end
		end
		n_it += 1
    end
	return lots
end

function VRP_iteratif(prp, q_dec, t, heuristic)
	n, k, Q = prp["n"], prp["k"], prp["Q"]
	cost, cc = calcul_dist(prp), [i for i in 1:n if q_dec[i,t] != 0]
	lots = switch_heuristic(prp, q_dec, t, heuristic)
	# TSP
	lots = TSP(lots, cost)
	
	lots = split_routing(lots, k)

	if length(lots) > k 
		lots = vehicles_splitting(lots, k, q_dec, t, Q, cost)
		if lots == -1
			return -1
		end
	end

	lots = lots_ordering(lots, cc, cost, q_dec, t, Q)

	lots = TSP2(lots, cost)

	filter!(e->e!=[0], lots) # RM EMPTY CYCLES

	return lots
end
