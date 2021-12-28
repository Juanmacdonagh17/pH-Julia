### A Pluto.jl notebook ###
# v0.17.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ dbd43667-6b86-4960-adfb-dbbd922a63ea
# begin
# 	using Optim, Plots, LaTeXStrings, PlutoUI, Unitful, PlotlyJS, DataFrames, StatsPlots
# 	plotly()
# end;
begin
	import Pkg
	    ENV["JULIA_MARGO_LOAD_PYPLOT"] = "no thank you"
	Pkg.activate(mktempdir())
	Pkg.add([
	    Pkg.PackageSpec(name="Plots", version="1"),
	    Pkg.PackageSpec(name="PlutoUI", version="0.7"),
	    Pkg.PackageSpec(name="Optim"),
	    Pkg.PackageSpec(name="LaTeXStrings"),
	    Pkg.PackageSpec(name="Unitful"),
	    Pkg.PackageSpec(name="DataFrames")
	    ])
	    
	using Plots
	using PlutoUI
	using LaTeXStrings
	using Optim
	using DataFrames
	plotly()
    end;

# ╔═╡ 5ad96e47-9ca9-493f-b899-a109826953b7
md"""
# pH solver en Julia

###### fuertisimamente basado en [pHsolve for python](https://github.com/rnelsonchem/pHcalc/blob/master/pHcalc/pHcalc.py)

"""



# ╔═╡ d136b0c6-5b94-11ec-01b6-4f5226fe2044
begin
mutable struct Acid
	ka::Union{Float64,Vector{Float64}}
	conc::Float64
	charge::Vector{Int64}
	pKA::Vector{Float64}
	_ka::Vector{Float64}
	function Acid(ka, conc, charge::Int64)
			if typeof(ka) <: Number
				ka=[ka]
			end

			_ka=copy(ka)
			sort!(ka, rev=:true)
			prepend!(_ka,1.)
			charge=charge:-1:charge-(length(ka)) |> collect
			# charge=LinRange(charge,-length(ka), length(ka))
			new(ka, conc,charge , -log10.(ka),_ka )
		end
end
	
mutable struct Neutral
	charge::Union{Vector{Int64},Int64}
	conc::Float64
end
	
struct System
	species
	function System(species...)
		new(species)
	end
end

function α(self::Acid,pH) # solo para un valor de pH, se podria pensar para un array de pHs.
		h3o=10.0^(-pH)
		power=length(self._ka)-1:-1:0
		h3o_power=(x->h3o^x).(power)
		Ka_prod=cumprod(self._ka)
		h3o_Ka=h3o_power.*(Ka_prod)
		h3o_Ka./sum(h3o_Ka)
	end

function α(self::Neutral,pH) # Devielve 1 porque no tiene disociacion, devuelvo en array para que se parezca a la otra implementación.
	[1]
end
	
function pHsolve(sys)
	function minimise(pH)
		h3o=10.0^(-pH)
     	oh = (10.0^(-14))/h3o
	 	x = (h3o - oh)
		for specie in sys.species
			x+=(specie.conc.*specie.charge.*α(specie,pH))|> sum
		end
	 	abs(x)
	end
	
	function phguess()
		phv=1:1.:14
		minimise.(phv) |> argmin |> x -> getindex(phv, x)
	end
	guess=phguess()
	optimize(x->minimise(first(x)), [guess], BFGS()).minimizer |> first
end
	T=Union{Acid, Neutral}
end
	

# ╔═╡ 38c1cace-25cc-4d66-830b-14d3c3aaba6d
begin
	fosforico=Acid([7.52e-3,6.23e-8,4.8e-13],.01,0)
	acetico=Acid(1.8e-5,.001,0)
	clorhidrico=Neutral(-1,1e-8)
	aspartico=Acid((x-> 10^-x).([2.09,9.82,3.86]),0.06, 1)
	sistema=System(clorhidrico)
end;

# ╔═╡ 1417b95f-d243-413e-a97f-38bc2af7fc32


# ╔═╡ 1e2f4ee6-0ec2-4357-adb5-8ab6c593d8df
especies=Dict(
	:aspartico=>aspartico,
	:fosforico=>fosforico,
	:acetico=>acetico,
	:clorhidrico=>clorhidrico
);

# ╔═╡ 184f4872-1992-40ff-b7c9-8ae35675a5e9
fosforico

# ╔═╡ f6d249c9-f32a-45ee-be3b-7a4a9764cfb8
begin
	phs=1:.1:14

	data=(x->α(fosforico,x)).(phs) |> x->hcat(x...) 
	# Plots.plot(phs,hcat(data...))
	p=Plots.plot()
	for i in 1:length(fosforico.charge)
		Plots.plot!(p, phs, data[i,:], label="especie $i")
		title!(p, "Diagrama de distribución de especies de H3PO4")
		xlabel!(p, "pH")
		ylabel!(p, "Fracción de Concentración")
	
		println(i)
	end
	for (i,pKa) in enumerate(fosforico.pKA)
		vline!(p, [pKa], label=:false)
end
	p
end

# ╔═╡ 635d433f-482d-45a9-a6fa-ab3b87b4b092

md"""

Sustancia a titular $(@bind esp Select(especies|> keys |> collect))


"""

# ╔═╡ 11691744-7892-4ac7-874b-e0c5dfb3c932
md"""
Concentración de titulante (M): $(@bind conc_na Slider(LinRange(.1,.5,50), default=0.1, show_value=true);)\

Volumen de ácido a titular (ml): $(@bind vol_erlen Slider(5:20, default = 10, show_value = true);)\

Concentración molar del Ácido (M): $(@bind conc_ac Slider(LinRange(.1,.5,50), show_value = true);)\


$(@bind tuki Button("Graficar!"))
"""

# ╔═╡ b29ab129-dbe6-40e1-bb86-53c6a7bf6c2c


# ╔═╡ 7b51c5c1-88e1-418b-851f-b728cc1a6a44
function pde(acido::T, conc_na, vol_erlen)
	times= typeof(acido) == Acid ? length(acido.ka) : 1 #Un acido fuerte solo puede perder un protón en ppio.
	vol = (chr -> (acido.conc * vol_erlen * chr)/conc_na)(collect(1:times))
	pH = map(1:times) do t
		base=Neutral(1, t* conc_na)
		System(acido, base) |> pHsolve
	end
	vol,pH
end;

# ╔═╡ 71afb541-c7de-40b8-bfa1-6a378ee54a6e
begin
	resolucion=60		
	# fosforico_tit=Acid([7.52e-3,6.23e-8,4.8e-13],conc_ac,0)
	sustancia=especies[esp]
	sustancia.conc=conc_ac
	vols_agregados=LinRange(.1,vol_erlen+5,resolucion)      #0.0001:.2:vol_erlen+5
	vols_en_erlen=vol_erlen .+ vols_agregados
	concs_en_erlen=(conc_na*cumsum(vols_agregados)) ./ vols_en_erlen
	pdes=pde(sustancia, conc_na, vol_erlen)
	pHsas=(x-> Neutral(1,x) |> x-> System(sustancia, x)|> pHsolve).(concs_en_erlen);
end;

# ╔═╡ 1708178e-c961-45d4-9dcb-5c763400c95a
begin
	tuki
	p2=Plots.scatter(vols_agregados ,pHsas, legend=:false)
	Plots.plot!(p2,vols_agregados ,pHsas, legend=:false )
	Plots.scatter!(pdes)
	xlabel!(p2,"Vol. NaOH")
	ylabel!("pH")
	title!("Titulación de Ácido $(esp) con NaOH")
	# Plots.savefig(p2,"falopa.png")
	p2;
	# png("faLoPa")
	# p2
end

# ╔═╡ 4d9fb42f-e7aa-4ffe-983f-4f8e90a1603f
pd=pde(fosforico, conc_na, vol_erlen)

# ╔═╡ 791fc74f-8535-44a4-831e-1a9b041619e2
Plots.scatter(pd)

# ╔═╡ aaad8512-d634-4cf1-82cd-923a276215fb
vol= (chr -> (conc_ac * vol_erlen * chr)/conc_na)([2,3])

# ╔═╡ d3eea87b-0695-4f39-a8b1-ef943619d094
@bind axa Select(especies |> keys |> collect)

# ╔═╡ f05cb0bf-12e1-48b0-9b3a-7d9f34018b6a
begin
	arr=1:.1:14
	p3 = (x-> α(especies[axa], x)).(arr) |> x->  hcat(x...) |> transpose |> x-> Plots.plot(arr,x, legend=:false)
	for pka ∈ especies[axa].pKA
		vline!([pka])
	end
	p3
end


# ╔═╡ 6fa6ead0-6512-4a5d-9d37-42876d15adc4
md"""
### Calculo vectorizado de pH *_(idea)_*
"""

# ╔═╡ 222595f9-5e92-48b3-862c-248c5f73ecfd
function minimise(sys::System,pH)
	h3o=10.0^(-pH)
     	oh = (10.0^(-14))/h3o
	 x = (h3o - oh)
	for specie in sys.species
		x+=(specie.conc.*specie.charge.*α(specie,pH))|> sum
	end
	 abs(x)
end

# ╔═╡ f29144ca-ebd5-403f-a4a8-70d136b43781
begin
	vec=0:.01:14
(x->minimise(sistema,x)).(vec) |> x-> Plots.plot(vec, x, yscale=:log) 
end

# ╔═╡ Cell order:
# ╟─5ad96e47-9ca9-493f-b899-a109826953b7
# ╠═dbd43667-6b86-4960-adfb-dbbd922a63ea
# ╠═d136b0c6-5b94-11ec-01b6-4f5226fe2044
# ╠═38c1cace-25cc-4d66-830b-14d3c3aaba6d
# ╠═1417b95f-d243-413e-a97f-38bc2af7fc32
# ╠═1e2f4ee6-0ec2-4357-adb5-8ab6c593d8df
# ╠═184f4872-1992-40ff-b7c9-8ae35675a5e9
# ╠═f6d249c9-f32a-45ee-be3b-7a4a9764cfb8
# ╟─635d433f-482d-45a9-a6fa-ab3b87b4b092
# ╠═11691744-7892-4ac7-874b-e0c5dfb3c932
# ╠═b29ab129-dbe6-40e1-bb86-53c6a7bf6c2c
# ╟─7b51c5c1-88e1-418b-851f-b728cc1a6a44
# ╟─71afb541-c7de-40b8-bfa1-6a378ee54a6e
# ╠═1708178e-c961-45d4-9dcb-5c763400c95a
# ╠═4d9fb42f-e7aa-4ffe-983f-4f8e90a1603f
# ╠═791fc74f-8535-44a4-831e-1a9b041619e2
# ╠═aaad8512-d634-4cf1-82cd-923a276215fb
# ╟─d3eea87b-0695-4f39-a8b1-ef943619d094
# ╠═f05cb0bf-12e1-48b0-9b3a-7d9f34018b6a
# ╟─6fa6ead0-6512-4a5d-9d37-42876d15adc4
# ╠═222595f9-5e92-48b3-862c-248c5f73ecfd
# ╠═f29144ca-ebd5-403f-a4a8-70d136b43781
