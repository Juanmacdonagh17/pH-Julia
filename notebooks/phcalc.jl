### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# ╔═╡ 5ad96e47-9ca9-493f-b899-a109826953b7
md"""
# pH solver en Julia

###### fuertisimamente basado en [pHsolve for python](https://github.com/rnelsonchem/pHcalc/blob/master/pHcalc/pHcalc.py)

"""



# ╔═╡ d136b0c6-5b94-11ec-01b6-4f5226fe2044
begin
struct Acid
	ka::Vector{Float64}
	charge::Vector{Float64}
	pKA::Vector{Float64}
	_ka::Vector{Float64}
	function Acid(ka)
			if typeof(ka) <: Number
				ka=[ka]
			end
			sort!(ka, rev=:true)
			_ka=copy(ka)
			append!(_ka,1.)
			charge=-collect(0:length(ka))
			new(ka, charge , -log10.(ka),_ka )
		end
end

function alpha(self::Acid,ph) # solo para un valor de pH, se podria pensar para un array de pHs
		h3o=10.0^(-ph)
		power=length(self._ka):-1:1
		h3o_power=(x->x^h3o).(power)
		Ka_prod=cumprod(self._ka)
		h3o_Ka=h3o_power.*(Ka_prod)
		h3o_Ka./sum(h3o_Ka)
	end
	
	
	
end


# ╔═╡ 38c1cace-25cc-4d66-830b-14d3c3aaba6d
begin
	fosforico=Acid([7.52e-3,6.23e-8,4.8e-13])
	alpha(fosforico, 7)## no funciona bien
end

# ╔═╡ Cell order:
# ╟─5ad96e47-9ca9-493f-b899-a109826953b7
# ╠═d136b0c6-5b94-11ec-01b6-4f5226fe2044
# ╠═38c1cace-25cc-4d66-830b-14d3c3aaba6d
