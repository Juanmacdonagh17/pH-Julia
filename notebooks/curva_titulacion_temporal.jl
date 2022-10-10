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

# ╔═╡ de7b7965-00e6-499d-8797-8588c4cda3ea
using Suppressor 

# ╔═╡ 4969064a-6f6d-11ec-1d9e-db3f5da90e7e
@suppress begin
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
		Pkg.PackageSpec(name="PHcalc"),
		Pkg.PackageSpec(name="OrderedCollections")
	    ])
	    
	using Plots, Statistics, PlutoUI, LaTeXStrings, Optim, DataFrames, PHcalc, OrderedCollections
	plotly()
    end;

# ╔═╡ 971751a5-b338-4e2a-80bd-c59e76ff3734
md"""

# Curvas de titulación


"""

# ╔═╡ 8fa44710-71a7-41f5-b1e2-878ad4368b7f
md"""
En esta notebook podemos explorar los comportamientos de las curvas de titulación para distintos ácidos, con concentraciones variables, frente a Hidróxido de Sodio, también con concentración variable. 


Sumado a esto, podemos analizar qué indicadores funcionan para determinar el punto de equivalencia de la misma, mediante un método gráfico.
"""

# ╔═╡ 4de35c9b-29bb-4e62-88d0-3e426df2aea5
begin
	md"""
	$\require{mhchem}$
	""";
	html"""<style>
main {
    max-width: 750px;
}
"""
end

# ╔═╡ ca743222-5da5-46d3-bfcd-53fab31acd8c
begin
	fosforico=Acid([7.52e-3,6.23e-8,4.8e-13],.01,0)
	acetico=Acid(1.8e-5,.001,0)
	clorhidrico=Neutral(-1,1e-8)
	aspartico=Acid((x-> 10^-x).([2.09,9.82,3.86]),0.06, 1)
	sistema=System(clorhidrico)
	especies=Dict(
		:Aspártico=>aspartico,
		:Fosfórico=>fosforico,
		:Acético=>acetico,
		:Clorhídrico=>clorhidrico
	);
	Indicadores=OrderedDict(
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
	),
		:Ninguno=>Nothing
)

	
	function pde(acido :: AbstractSpecies, conc_na, vol_ac)
		times= typeof(acido) == Acid ? length(acido.ka) : 1 #Un acido fuerte solo puede perder un protón en ppio.
		vol = (chr -> (acido.conc * vol_ac * chr)/conc_na)(collect(1:times))
		conc_ac=acido.conc
		# Chr representa las distintas cargas que puede tener un ácido poliprótico, por ejemplo para Ácido fosfórico, las concentraciones de los puntos de equivalencia deberían ser obtenerse multiplicando [1,2,3]*Concentración_base*volumen_base
		pH = map( vol ) do v
			acido.conc = ( conc_ac * vol_ac ) / ( v + vol_ac ) # a)
			base=Neutral(+1, #carga
				 ( ( conc_na * v )/( v + vol_ac ) )  #b)
			)
			System(acido, base) |> s -> pHfast(s,.1)
		end
		
		(
			vol=vol,
		 	pH=pH
		)

	end
end;


#en a) y b) corrijo por dilución 

# ╔═╡ 5052e33a-687b-43e4-8437-b8fe392cea87
md"""

## Explicación teórica:

Para entender las curvas de pH podemos tratarlo como una funcion del volumen de Base:

$$\text{pH} (\text{Vol}_{\ce{NaOH}})$$

Si recordemos la definición de pH: 

$$\text{pH}=-  \log{([H^+])}$$

Podemos observar que, por cada gota de base que agregamos, modificamos la concentracion de protones por la siguiente reaccion: 

$$\ce{H+ +OH- <=>H2O}$$

Entonces, a tiempo $${0}$$ (cuando no hay ninguna gota de base), el pH inicial se calcula como: 

$$\text{pH}_{V_0} = -\log{[H^+_i]}$$


Pero al agregar la primer gota de Base, este se modifica de la siguiente manera:


$$\text{pH}_{V_1}=-\log{\Bigg(\dfrac{(\text{moles }H^{+}_{0}- \text{ moles }{OH}^-_{1})*(V_i+V_{1\text{NaOH}})}{1000\text{ml}}\Bigg)}$$


En esta ecuacion le restamos a los moles *iniciales* (en el punto anterior) debido a que modifican las concentraciones.

De forma general, podemos expresar la función de la siguiente manera:

$$\text{pH}_{V_x}=-\log{\Bigg(\dfrac{(\text{moles }H^{+}_{x-1}-\text{moles }{OH}^-_{1})*(V_{x-1}+V_{x(\text{NaOH})})}{1000\text{ml}}\Bigg)}$$

Con esto ya tenemos (casi) todo lo necesario para armar nuestra función. Los parametros a usar van a ser:\

1) Concentración inicial del titulante
2) Volumen de ácido 
3) Concentración inicial de Ácido


"""


# ╔═╡ a181e191-ec53-49f3-b122-87aa80a73074
md"""
## Parámetros a graficar:
"""

# ╔═╡ f56b3295-d66c-4fbd-ba83-4c54e9d6e871
md"""


Ácido a titular $(@bind esp Select(especies|> keys |> collect))

Concentración de titulante (M): $(@bind conc_na Slider(LinRange(.1,.5,50), default=0.1, show_value=true);)\

Volumen de ácido a titular (ml): $(@bind vol_ac Slider(5:20, default = 10, show_value = true);)\

Concentración molar del Ácido (M): $(@bind conc_ac Slider(LinRange(.1,.5,50), show_value = true);)\
\

Indicador: $(@bind ind Select(collect(keys(Indicadores))))

"""

# ╔═╡ 36d23bf8-99ce-41f9-84d0-d59da6ab44b5
begin
	especie=especies[esp] 
	resolucion=200
	vol_na= 3*(conc_ac*vol_ac)/conc_na
	vol_agregados=LinRange(0,vol_na + 10, resolucion)
	conc_na_erlen=(conc_na*vol_agregados)./(vol_agregados.+vol_ac)
	conc_ac_erlen=(conc_ac*vol_ac)./(vol_agregados.+vol_ac)
	especie.conc=conc_ac
	punto_de_equivalencia=pde(especie,conc_na ,vol_ac)
	
	
	pH=map(conc_na_erlen, conc_ac_erlen) do Cₙ,Cₜ
		especie.conc=Cₜ
		System(
			Neutral(1, Cₙ), 
			especie
		) |> s -> pHfast(s,.1)
	end
end;

# ╔═╡ 371acdab-ee5e-4adc-b373-87c3fe78e8af
md"""
## Gráfico:
"""

# ╔═╡ e9ae7938-e4cf-4a28-932e-5b956436b53d
begin
	p = scatter(vol_agregados,pH, legend=:false, size= (1000,600))
	plot!(vol_agregados, pH, label=:false)
	scatter!(punto_de_equivalencia.vol, punto_de_equivalencia.pH, markershape=:hexagon, markersize=12 ,markercolor=:green)
	indicador=Indicadores[ind]
	if ind ≠ :Ninguno
		hline!([indicador[:viraje][1]], color = indicador[:colores][1], width = 2)
		hline!([indicador[:viraje][2]], color = indicador[:colores][2], width = 2)
	end
	xlabel!("Volumen de NaOH (ml)")
	ylabel!("pH")
end

# ╔═╡ fe6e663b-3444-4b0f-9b1f-b3694aaec3e5
TableOfContents(title="Curva de titulación Ácido Base 🎢")

# ╔═╡ 4cd4972b-6886-4217-a4a8-42527f253e49
md"""
Desarrollado por Eduardo Gonik y Juan Mac Donagh. UNLP, Facultad de Ciencias Exactas, 2022.
"""

# ╔═╡ Cell order:
# ╟─971751a5-b338-4e2a-80bd-c59e76ff3734
# ╟─8fa44710-71a7-41f5-b1e2-878ad4368b7f
# ╟─de7b7965-00e6-499d-8797-8588c4cda3ea
# ╟─4969064a-6f6d-11ec-1d9e-db3f5da90e7e
# ╟─4de35c9b-29bb-4e62-88d0-3e426df2aea5
# ╟─ca743222-5da5-46d3-bfcd-53fab31acd8c
# ╟─5052e33a-687b-43e4-8437-b8fe392cea87
# ╟─a181e191-ec53-49f3-b122-87aa80a73074
# ╟─f56b3295-d66c-4fbd-ba83-4c54e9d6e871
# ╟─36d23bf8-99ce-41f9-84d0-d59da6ab44b5
# ╟─371acdab-ee5e-4adc-b373-87c3fe78e8af
# ╟─e9ae7938-e4cf-4a28-932e-5b956436b53d
# ╟─fe6e663b-3444-4b0f-9b1f-b3694aaec3e5
# ╟─4cd4972b-6886-4217-a4a8-42527f253e49
