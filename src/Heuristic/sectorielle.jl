function SECT(data, q_dec, t) 
	Q = data["Q"]
	n, Q = data["n"], data["Q"]
	coord = [(data["Clients"].X[i], data["Clients"].Y[i]) for i in 1:data["n"]+1] # Clean coords to get easy manipulated tuples /!\ we include the factory
    mid, r = [coord[1][1], coord[1][2]], [1, 0]
    cc, c_ang = [], []

	for i in 1:n
		if q_dec[i,t] != 0
			push!(cc, i)
		end
	end
	for i in cc
		pt_mid = [coord[i+1][1], coord[i+1][2]] - mid
		deg = mod(rad2deg(atan(r...) - atan(pt_mid...)), 360)
		push!(c_ang, (i, deg))
	end
	sort!(c_ang, by = x->x[2]) 
	cycles, cycle, _sum = Vector{Int64}[], [0], 0
	for ang in c_ang
		c = ang[1]
		_sum += q_dec[c, t]
		if _sum > Q
			push!(cycles, cycle)
			_sum, cycle = q_dec[c, t], [0]
		end
		push!(cycle, c)
	end
	if length(cycle) > 1
		push!(cycles, cycle)
	end
	return cycles
end