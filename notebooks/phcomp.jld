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

# ╔═╡ fe74efb7-5fe3-4bb4-896a-a8ebc8c5a74e

begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
	    Pkg.PackageSpec(name="Plots", version="1"),
	    Pkg.PackageSpec(name="PlutoUI", version="0.7"),
	    Pkg.PackageSpec(name="Optim"),
	    Pkg.PackageSpec(name="LaTeXStrings"),
	    Pkg.PackageSpec(name="Unitful"),
	    Pkg.PackageSpec(name="DataFrames"),
		 Pkg.PackageSpec(name="Statistics")
	    ])
	    
	using Plots
	using Statistics
	using PlutoUI
	using LaTeXStrings
	using Optim
	using DataFrames
	plotly()
    end;

# ╔═╡ 77dec32e-6a8c-11ec-0c99-bfa1f7fe5392
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


function pHfast(sys, precision=.01)
	function minimise(pH)
		h3o=10.0^(-pH)
     	oh = (10.0^(-14))/h3o
	 	x = (h3o - oh)
		for specie in sys.species
			x+=(specie.conc.*specie.charge.*α(specie,pH))|> sum
		end
	 	abs(x)
	end
	phv=1:precision:14
	minimise.(phv) |> argmin |> x -> getindex(phv,x) 
		
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
	

# ╔═╡ 1bd376ba-c267-4585-92cd-dff4a3a6658e
md"""
# Comparación de métodos pHcalc y pHfast (con precision .1 y .01)


"""

# ╔═╡ f4ac5988-df2d-47eb-a798-031490cf837f
md"""
El enfoque de este notebook es comparar los distintos métodos de calcular pH usando la lib pHcalc.jl. La comparacion que se muestra es solamente hecha en cuanto a la precisión, no velocidad. Esto último es porque en Pluto las celdas se ejecutan secuencialmente cuando ocurren cambios que las afectan, entonces no tiene sentido comparar los plots con el mismo set de sliders.

"""

# ╔═╡ c8c17787-3452-41a2-8928-2d37f2016243
begin
	fosforico=Acid([7.52e-3,6.23e-8,4.8e-13],.01,0)
	acetico=Acid(1.8e-5,.001,0)
	clorhidrico=Neutral(-1,1e-8)
	aspartico=Acid((x-> 10^-x).([2.09,9.82,3.86]),0.06, 1)
	sistema=System(clorhidrico)
end;

# ╔═╡ f68efe0b-2a45-400f-81be-dfb6c5aabc48
especies=Dict(
	:aspartico=>aspartico,
	:fosforico=>fosforico,
	:acetico=>acetico,
	:clorhidrico=>clorhidrico
);

# ╔═╡ 9ea0cec9-acaa-4c12-a0af-5e352304e79f
function pde(acido::T, conc_na, vol_erlen)
	times= typeof(acido) == Acid ? length(acido.ka) : 1 #Un acido fuerte solo puede perder un protón en ppio.
	vol = (chr -> (acido.conc * vol_erlen * chr)/conc_na)(collect(1:times))
	pH = map(1:times) do t
		base=Neutral(1, t* conc_na)
		System(acido, base) |> pHsolve
	end
	vol,pH
end;

# ╔═╡ d3e4683a-9d8e-4375-b621-089863dac270

md"""

Sustancia a titular $(@bind esp Select(especies|> keys |> collect))

Concentración de titulante (M): $(@bind conc_na Slider(LinRange(.1,.5,50), default=0.1, show_value=true);)\

Volumen de ácido a titular (ml): $(@bind vol_erlen Slider(5:20, default = 10, show_value = true);)\

Concentración molar del Ácido (M): $(@bind conc_ac Slider(LinRange(.1,.5,50), show_value = true);)\


"""

# ╔═╡ 7d468210-f0e5-482b-8891-d993582b24b7
begin
	resolucion=60		
	# fosforico_tit=Acid([7.52e-3,6.23e-8,4.8e-13],conc_ac,0)
	sustancia=especies[esp]
	sustancia.conc=conc_ac
	vols_agregados=LinRange(.1,vol_erlen+5,resolucion)      #0.0001:.2:vol_erlen+5
	vols_en_erlen=vol_erlen .+ vols_agregados
	concs_en_erlen=(conc_na*cumsum(vols_agregados)) ./ vols_en_erlen
	pdes=pde(sustancia, conc_na, vol_erlen)
	phfast01=(x-> Neutral(1,x) |> x-> System(sustancia, x)|> x-> pHfast(x,.1)).(concs_en_erlen);
	phfast001=	phfast01=(x-> Neutral(1,x) |> x-> System(sustancia, x)|> x-> pHfast(x,.01)).(concs_en_erlen);
	phcalc=	phfast01=(x-> Neutral(1,x) |> x-> System(sustancia, x)|> pHsolve).(concs_en_erlen);
end;

# ╔═╡ d34796b1-bd57-41e5-98ca-9d7de782ce2b
begin
	p=Plots.scatter(vols_agregados ,phfast01,label="fast 0.1")
	Plots.scatter!(vols_agregados ,phfast001,label="fast 0.01")
	Plots.scatter!(vols_agregados, phcalc, label="metodo completo")
	Plots.plot!(vols_agregados ,phfast01, label=:false)
	Plots.plot!(vols_agregados ,phfast001, label=:false)
	Plots.plot!(vols_agregados ,phcalc, label=:false)
	
	Plots.scatter!(pdes, label=:false)
	xlabel!("Vol. NaOH")
	ylabel!("pH")
	title!("Titulación de Ácido $(esp) con NaOH")
	p;
end

# ╔═╡ Cell order:
# ╟─fe74efb7-5fe3-4bb4-896a-a8ebc8c5a74e
# ╟─77dec32e-6a8c-11ec-0c99-bfa1f7fe5392
# ╟─1bd376ba-c267-4585-92cd-dff4a3a6658e
# ╟─f4ac5988-df2d-47eb-a798-031490cf837f
# ╟─c8c17787-3452-41a2-8928-2d37f2016243
# ╟─f68efe0b-2a45-400f-81be-dfb6c5aabc48
# ╟─9ea0cec9-acaa-4c12-a0af-5e352304e79f
# ╟─d3e4683a-9d8e-4375-b621-089863dac270
# ╟─7d468210-f0e5-482b-8891-d993582b24b7
# ╟─d34796b1-bd57-41e5-98ca-9d7de782ce2b
