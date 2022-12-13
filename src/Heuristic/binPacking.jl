function BP(data, q_dec, t)
    n, Q, _sum, lots, cur = data["n"], data["Q"], 0, Vector{Int64}[], [0]
	for i in 1:n 
		if q_dec[i,t] != 0
			_sum += q_dec[i,t]
			if _sum > Q
				push!(lots, cur)
				_sum = q_dec[i,t]
				cur = [0]
			end
			push!(cur, i)
		end
	end
	if length(cur) > 1
		push!(lots, cur)
	end
	return lots
end