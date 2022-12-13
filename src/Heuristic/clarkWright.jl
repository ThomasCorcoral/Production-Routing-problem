include("../utils.jl")

function CW(prp, q_dec, t)
	# cc -> cost customer 
	n, cc, d, res, s = prp["n"], [], 0, Array[], []
	cost = calcul_dist(prp)
	for i in 1:n
		if q_dec[i,t] != 0
			d += q_dec[i,t]
			push!(cc, i)
		end
	end
	for i in cc
		push!(res, [0,i])
	end
	for i in cc 
		for j in cc 
			if i!=j
				push!(s, [(i, j), cost[1, i+1] + cost[1, j+1] - cost[i+1, j+1] + abs(cost[1, i+1] + cost[1, j+1]) + (q_dec[i, t]+q_dec[j, t]) / d])
			end
		end
	end
	sort!(s, by = x->x[2], rev=true) 
	for k in 1:length(s)
		i, j = s[k][1][1], s[k][1][2]
		b_i, b_j, i_b_i, i_b_j = [], [], -1, -1
		for ii in 1:length(res) 
			if j in res[ii]
				b_i, i_b_i = res[ii], ii
			end
			if j in res[ii]
				b_j, i_b_j = res[ii], ii
			end
			if i_b_i != -1 && i_b_j != -1
				break
			end
		end
		if i_b_j != i_b_i 
			fused = copy(b_i) 
			for e in b_j
				if !(e in fused)
					push!(fused, e)
				end
			end
			fuse = 0 
			for l in 2:length(fused)
				fuse += q_dec[fused[l], t]
			end
			if fuse <= prp["Q"]
				new_res = []
				for ii in 1:length(res)
					if ii == i_b_i || ii == i_b_j 
						continue 
					end
					push!(new_res, res[ii])
				end
				res = new_res
				push!(res, fused)
			end
		end
	end
	return res 
end
