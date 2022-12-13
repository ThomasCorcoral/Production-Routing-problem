using GraphPlot
using Cairo
using Compose
using Fontconfig
using Colors

include("utils.jl")

function graph_to_png(G, clients_t, qt, e_col, n_col, coordx_t, coordy_t, filename)
	gp = gplot(G, 
			coordx_t, 
			coordy_t, 
			nodelabel = [(clients_t[i], qt[i]) for i in 1:length(clients_t)], 
			edgelinewidth = [200*length(clients_t) for e in edges(G)], 
			edgestrokec = e_col, 
			nodefillc = n_col)
	draw(PNG(filename, 3*length(clients_t)cm, 3*length(clients_t)cm), gp) 
end

function generate_colors(G, clients_t, circuits, d)
	n_col = [colorant"red" for i in 1:length(clients_t)]
	dictEdge = Dict()
	nc = 1
	cols = distinguishable_colors(length(circuits), [RGB(0,0,0.8), RGB(0.8,0,0)])
	for circuit in circuits
		n_col[1] = RGBA(1,1,1) # The factory is WHITE
		for (i, j) in zip(circuit[1:end-1], circuit[2:end])
			# println(i, "/", j)
			# println(d)
			# println(d[i], "/", d[j])
			# println()
			add_edge!(G, d[i], d[j]) 
			dictEdge[(d[i], d[j])] = cols[nc]
			n_col[d[j]] = cols[nc]
		end

		add_edge!(G, d[circuit[end]], 1) 
		dictEdge[(d[circuit[end]], 1)] = cols[nc]
		nc += 1
	end
	e_col = []
	for (e_idx, e) in enumerate(edges(G))
        i = src(e)
        j = dst(e)
		push!(e_col, dictEdge[(i,j)])
	end
	return G, e_col, n_col
end

function pdi_to_png(prp, q, t, circuits, filename, date)

	# println(prp)
	# println()
	# println(q)
	# println()
	# println(t)
	# println()
	# println(circuits)
	# println()

	clients_t = [0] # Init for factory
	for i in 1:prp["n"]
		if q[i, t] != 0
			push!(clients_t, i)
		end
	end

	d = Dict() # element : index
	for ind in 1:length(clients_t)
		if !(clients_t[ind] in keys(d))
			d[clients_t[ind]] = ind
		end
	end

	G = DiGraph(length(clients_t))
	coordx_t, coordy_t = generate_coords(prp, d, clients_t)

	G, e_col, n_col = generate_colors(G, clients_t, circuits, d)
	
	qt = [0.0] # init for factory
	for i in clients_t[2:length(clients_t)]
		push!(qt, floor.(q[i,t]))
	end

	if !isdir("log")
		mkdir("log")
	end
	path = string("log/", filename)
	if !isdir(path)
		mkdir(path)
	end
	graph_to_png(G, clients_t, qt, e_col, n_col, coordx_t, coordy_t, path * "/" * "t_" * string(t) * ".png")
end