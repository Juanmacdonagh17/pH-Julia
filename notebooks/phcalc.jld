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
		 Pkg.PackageSpec(name="Statistics"),
		Pkg.PackageSpec(url="https://github.com/egonik-unlp/PHcalc.jl")
	    ])
	    
	using Plots, Statistics, PlutoUI, LaTeXStrings, Optim, DataFrames, PHcalc
	plotly()
    end;

# ╔═╡ 5ad96e47-9ca9-493f-b899-a109826953b7
md"""
# pH solver en Julia

###### fuertisimamente basado en [pHsolve for python](https://github.com/rnelsonchem/pHcalc/blob/master/pHcalc/pHcalc.py)

"""



# ╔═╡ 38c1cace-25cc-4d66-830b-14d3c3aaba6d
begin
	fosforico=Acid([7.52e-3,6.23e-8,4.8e-13],.01,0)
	acetico=Acid(1.8e-5,.001,0)
	clorhidrico=Neutral(-1,1e-8)
	aspartico=Acid((x-> 10^-x).([2.09,9.82,3.86]),0.06, 1)
	sistema=System(clorhidrico)
end;

# ╔═╡ 1e2f4ee6-0ec2-4357-adb5-8ab6c593d8df
especies=Dict(
	:aspartico=>aspartico,
	:fosforico=>fosforico,
	:acetico=>acetico,
	:clorhidrico=>clorhidrico
);

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

# ╔═╡ 7b51c5c1-88e1-418b-851f-b728cc1a6a44
function pde(acido::T, conc_na, vol_erlen) where T<: Union{Acid, Neutral}
	times= typeof(acido) == Acid ? length(acido.ka) : 1 #Un acido fuerte solo puede perder un protón en ppio.
	vol = (chr -> (acido.conc * vol_erlen * chr)/conc_na)(collect(1:times))
	pH = map(1:times) do t
		base=Neutral(1, t* conc_na)
		System(acido, base) |> pHsolve
	end
	vol,pH
end;

# ╔═╡ 8e603414-5e4f-42c6-bb0c-63538016c50d
begin
	resolucion=60		
	# fosforico_tit=Acid([7.52e-3,6.23e-8,4.8e-13],conc_ac,0)
	sustancia=especies[esp]
	sustancia.conc=conc_ac
	vols_agregados=LinRange(.1,vol_erlen+5,resolucion)      #0.0001:.2:vol_erlen+5
	vols_en_erlen=vol_erlen .+ vols_agregados
	concs_en_erlen=(conc_na*vols_agregados) ./ vols_en_erlen
	pdes=pde(sustancia, conc_na, vol_erlen)
	pHsas=(x-> Neutral(1,x) |> x-> System(sustancia, x)|> x-> pHfast(x,.1)).(concs_en_erlen);
end;

# ╔═╡ 4d9fb42f-e7aa-4ffe-983f-4f8e90a1603f
begin
	tuki
	p5=Plots.scatter(vols_agregados ,pHsas,label=:false)
	Plots.plot!(vols_agregados ,pHsas, label="fast" )
	Plots.scatter!(pdes, label=:false)
	xlabel!("Vol. NaOH")
	ylabel!("pH")
	title!("Titulación de Ácido $(esp) con NaOH")
	p5;
end

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


# ╔═╡ Cell order:
# ╟─5ad96e47-9ca9-493f-b899-a109826953b7
# ╠═dbd43667-6b86-4960-adfb-dbbd922a63ea
# ╠═38c1cace-25cc-4d66-830b-14d3c3aaba6d
# ╠═1e2f4ee6-0ec2-4357-adb5-8ab6c593d8df
# ╟─f6d249c9-f32a-45ee-be3b-7a4a9764cfb8
# ╠═635d433f-482d-45a9-a6fa-ab3b87b4b092
# ╠═11691744-7892-4ac7-874b-e0c5dfb3c932
# ╠═7b51c5c1-88e1-418b-851f-b728cc1a6a44
# ╠═8e603414-5e4f-42c6-bb0c-63538016c50d
# ╟─4d9fb42f-e7aa-4ffe-983f-4f8e90a1603f
# ╟─d3eea87b-0695-4f39-a8b1-ef943619d094
# ╟─f05cb0bf-12e1-48b0-9b3a-7d9f34018b6a
