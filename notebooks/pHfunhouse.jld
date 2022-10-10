### A Pluto.jl notebook ###
# v0.19.12

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

# ╔═╡ 60bc58bb-7c95-484b-a69e-3de41f022f97
begin
	import Pkg
	    ENV["JULIA_MARGO_LOAD_PYPLOT"] = "no thank you"
	Pkg.activate(mktempdir())
	Pkg.add([
	    Pkg.PackageSpec(name="Plots", version="1"),
	    Pkg.PackageSpec(name="PlutoUI", version="0.7"),

	    Pkg.PackageSpec(name="LaTeXStrings")


	    ])
	    
	using Plots
	using PlutoUI
	using LaTeXStrings
	
	plotly()
    end;


# ╔═╡ a9aa8b20-48a6-11ec-037f-e1686702dd38
md"""
# Curva de titulación Ácido-Base en Julia: 💻
"""

# ╔═╡ 78375034-5f6b-4055-a373-fe40106764b9
md"""
$\require{mhchem}$
$\require{relsize}$
""";

# ╔═╡ 18d6f3bc-f18e-4e40-9aa4-92d25bb467f6


# ╔═╡ 4ec4df4b-b44d-474c-92ef-4e4cd1cd375f
#PlutoSliderServer.run_notebook(path_to_notebook)

# ╔═╡ 81034c90-049c-44b9-967b-23a5ce14ac98
md"""
### Explicación teórica (para ácido y base *fuertes*)
"""

# ╔═╡ 1fe00159-7a68-4747-a2b0-8ae4fb3d911f
md"""
Para graficas hay que tratar al pH como una funcion del volumen de Base:

$$\text{pH} (\text{Vol}_{\ce{NaOH}})$$

Recordemos como se define el pH: 

$$\text{pH}=-  \log{([H^+])}$$

Por cada gota de base que agregamos, modificamos la concentracion de protones por la siguiente reaccion: 

$$\ce{H+ +OH- <=>H2O}$$

Entonces, en el momento 0, cuando no hay ninguna gota de base, el pH inicial se calcula como: 

$$\text{pH}_{V_0} = -\log{[H^+_i]}$$


Pero al agregar la primer gota de Base, este se modifica de la siguiente manera:


$$\text{pH}_{V_1}=-\log{\Bigg(\dfrac{(\text{moles }H^{+}_{0}- \text{ moles }{OH}^-_{1})*(V_i+V_{1\text{NaOH}})}{1000\text{ml}}\Bigg)}$$


En esta ecuacion le restamos a los moles *iniciales* (en el punto anterior) debido a que modifican las concentraciones.

De forma general, podemos expresar la función de la siguiente manera:

$$\text{pH}_{V_x}=-\log{\Bigg(\dfrac{(\text{moles }H^{+}_{x-1}-\text{moles }{OH}^-_{1})*(V_{x-1}+V_{x(\text{NaOH})})}{1000\text{ml}}\Bigg)}$$

Con esto ya tenemos (casi) todo lo necesario para armar nuestra función. Los parametros a usar van a ser (para titular ácidos y bases fuertes!):\

1) Concentración inicial de Ácido
2) Volumen inicial de ácido 
3) Concentración inicial de Base
4) Valores de base agregados
"""

# ╔═╡ 6fac53cd-9441-4c65-9c4b-a24c10449aa4
md"""
#
"""

# ╔═╡ af14eab0-655e-425b-a483-dc63fc501045
md"""
#
"""

# ╔═╡ c05192ea-46d5-4a4d-8d7b-aabd5035c241
Indicadores=Dict(
	"Fenolftaleína"=>
Dict(
		:viraje=> [10,8.2],
		:colores=>[:pink,:whitesmoke]
),
	"Naranja de metilo"=>
	Dict(
		:viraje=>[4.4,3.1],
		:colores=>[:red,:yellow]
	),
	"Alizarina"=>
	Dict(
		:viraje=>[12.4,11],
		:colores=>[:red,:yellow]
	),
	"Rojo de cresol"=>
	Dict(
		:viraje=>[8.8,7],
		:colores=>[:yellow,:red]
	),
	"Violeta de metilo"=>
	Dict(
		:viraje=>[1.6,.2],
		:colores=>[:blueviolet,:yellow]
		)
);


# ╔═╡ 8b512d53-1a96-40f8-88d6-9caf28bbe9ef
md"""
Sustancias a usar (esto no hace nada todavía):\
$(@bind e Select(["first" => "HCL", "second"=>"HAc"])) y
$(@bind d Select(["first" => "NaOH", "second"=>"Base debil"]))
"""

# ╔═╡ 1aefcb86-c78f-4bee-9205-75c2b1d972d3
md"""
Elegimos los valores para nuestra titulación: \
Volumen total de la bureta (ml): $(@bind w Slider(11:50, default = 25, show_value = true))\
Volumen de ácido a titular (ml): $(@bind x Slider(5:20, default = 10, show_value = true))\
Concentración molar del Ácido: $(@bind y Slider(0.1: 0.05 :1, show_value = true))\
Concentración molar de la Base: $(@bind z Slider(0.1: 0.05 :1, show_value = true))\
"""

# ╔═╡ e803cbfc-b47c-4168-9ca7-19656e7f7202
begin
	bureta = collect(1:0.2:w)
	erlen = x
	erlen_t = []
	for i in range(1,length=length(bureta))
		if i == 1
			append!(erlen_t,erlen+(bureta[i]))
		else
			append!(erlen_t,(erlen+bureta[i]))
		end
	end
end

# ╔═╡ 8f64a5cd-c290-4972-8714-f914d730f700
begin
	mol_HCL = (y*x)/1000
	con_NaOH = z
	pH_list = []
	con_list = []
	for i in range(1,length=length(bureta))
		append!(con_list,(((mol_HCL-((con_NaOH*bureta[i])/1000))*(1000))/erlen_t[i]))
		if con_list[i] > 0
			append!(pH_list,(-log10(con_list[i])))
		else  
			append!(pH_list,(14+log10(-con_list[i])))
		end
	end
	#for i in range(1,length=length(pH_list))
		#if pH_list[i] == -Inf
			#pH_list_lie = replace(pH_list,pH_list[i] =>(pH_list[i-1]+pH_list[i+1])/2)
			#pH_list[i] == (pH_list[i-1]+pH_list[i+1])/2
		#end
	#end
end

# ╔═╡ 445a3f83-f82b-4175-a85c-f84c3bb09d3a
begin 
	scatter(bureta,pH_list, title = "Curva de pH", label = ["pH" "pH"], xlabel = "ml NaOH", ylabel = "pH", legend = false, mode="lines");
	if last(pH_list) > 7
		if -Inf in pH_list 
			vol_eq = bureta[findall(isequal(-Inf), pH_list)]
			Plots.plot!((bureta,pH_list), title= "Curva de pH", legend =false, label = false, hover=false);
			Plots.plot!((vol_eq,7), title= "Curva de pH", legend =false, markercolor = :red, markershape = :pentagon, label = "pH 7");
		else
			v1 = bureta[findall(isequal(last((pH_list[pH_list .< 3]))),pH_list)]
			v2 = bureta[findall(isequal(first((pH_list[pH_list .> 3]))),pH_list)]
			vol_eq = (v1.+v2)./2
			Plots.plot!((bureta,pH_list), title= "Curva de pH", legend =false, label = false, hover=false);
			Plots.plot!((vol_eq,7), title= "Curva de pH", legend =false, markercolor = :red, markershape = :pentagon, label = "pH 7");
		end
	else
		Plots.plot!((bureta,pH_list), title= "Curva de pH", legend =false)
	end
end

# ╔═╡ a6f6bacb-f8c3-4bca-9b85-964afb90f662
md"""
Indicadores:
"""

# ╔═╡ 537ea4b1-8cd1-4897-97be-a859d6519b68
@bind ind Select(collect(keys(Indicadores)))

# ╔═╡ e7a769e6-7cbd-4aa3-af4f-eea1dbb34d45
begin
	indicador=Indicadores[ind];
	p2 = scatter(bureta,pH_list, title = "Curva de pH", label = ["pH" "pH"], xlabel = "ml NaOH", ylabel = "pH", legend = false);
	
	hline!(p2,[indicador[:viraje][1]], color = indicador[:colores][1], width = 2);
	hline!(p2,[indicador[:viraje][2]], color = indicador[:colores][2], width = 2)
end

# ╔═╡ ebf4b542-dcd8-41c9-b567-6433abef0e87


# ╔═╡ 3d7010d3-3ced-4fe6-bcce-01e41d656c57
#if pH_list[46] == -Inf
	#a = replace(pH_list,pH_list[46] =>(pH_list[46-1]+pH_list[46+1])/2)
			#pH_list[i] == (pH_list[i-1]+pH_list[i+1])/2
#end

# ╔═╡ 6ff68f15-b5fa-490c-84c0-115997008326
md"""
Number of  HO's $(@bind HO Slider(1:1:10))

Number of  H's $(@bind H Slider(1:1:10))

"""  

# ╔═╡ d3fb6980-0ab9-4b98-b7fc-3928a733af91
L"""
\ce{%$H H+ + %$HO HO- <=> H2O } 
"""

# ╔═╡ Cell order:
# ╟─a9aa8b20-48a6-11ec-037f-e1686702dd38
# ╟─78375034-5f6b-4055-a373-fe40106764b9
# ╠═18d6f3bc-f18e-4e40-9aa4-92d25bb467f6
# ╟─60bc58bb-7c95-484b-a69e-3de41f022f97
# ╟─4ec4df4b-b44d-474c-92ef-4e4cd1cd375f
# ╟─81034c90-049c-44b9-967b-23a5ce14ac98
# ╠═1fe00159-7a68-4747-a2b0-8ae4fb3d911f
# ╟─6fac53cd-9441-4c65-9c4b-a24c10449aa4
# ╠═af14eab0-655e-425b-a483-dc63fc501045
# ╠═c05192ea-46d5-4a4d-8d7b-aabd5035c241
# ╟─8b512d53-1a96-40f8-88d6-9caf28bbe9ef
# ╠═1aefcb86-c78f-4bee-9205-75c2b1d972d3
# ╠═e803cbfc-b47c-4168-9ca7-19656e7f7202
# ╠═8f64a5cd-c290-4972-8714-f914d730f700
# ╠═445a3f83-f82b-4175-a85c-f84c3bb09d3a
# ╟─a6f6bacb-f8c3-4bca-9b85-964afb90f662
# ╠═537ea4b1-8cd1-4897-97be-a859d6519b68
# ╠═e7a769e6-7cbd-4aa3-af4f-eea1dbb34d45
# ╠═ebf4b542-dcd8-41c9-b567-6433abef0e87
# ╟─3d7010d3-3ced-4fe6-bcce-01e41d656c57
# ╟─6ff68f15-b5fa-490c-84c0-115997008326
# ╠═d3fb6980-0ab9-4b98-b7fc-3928a733af91
