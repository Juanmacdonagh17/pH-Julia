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

# ╔═╡ 4969064a-6f6d-11ec-1d9e-db3f5da90e7e

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
		Pkg.PackageSpec(name="PHcalc"),
		Pkg.PackageSpec(name="OrderedCollections")
	    ])
	    
	using Plots, Statistics, PlutoUI, LaTeXStrings, Optim, DataFrames, PHcalc, OrderedCollections
	plotly()
    end;

# ╔═╡ 4de35c9b-29bb-4e62-88d0-3e426df2aea5
begin
	md"""
	$\require{mhchem}$
	""";
	html"""<style>
main {
    max-width: 1500px;
}
"""
end

# ╔═╡ 971751a5-b338-4e2a-80bd-c59e76ff3734
md"""

# Curva de titulación



Cree esta notebook explicando el código de la curva de titulación que no me doy cuenta que problema tiene, pero me parece que está dando mal. También está el código para calcular los puntos de equivalencia, que todavía no me fijé si está andando bien

"""

# ╔═╡ ca743222-5da5-46d3-bfcd-53fab31acd8c
begin
	fosforico=Acid([7.52e-3,6.23e-8,4.8e-13],.01,0)
	acetico=Acid(1.8e-5,.001,0)
	clorhidrico=Neutral(-1,1e-8)
	aspartico=Acid((x-> 10^-x).([2.09,9.82,3.86]),0.06, 1)
	sistema=System(clorhidrico)
	especies=Dict(
		:aspartico=>aspartico,
		:fosforico=>fosforico,
		:acetico=>acetico,
		:clorhidrico=>clorhidrico
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

### Curva de titulación

Para la curva de titulación partimos de los siguientes parámetros:


$V_{acido}, C_{acido}, C_{base}$

Además de $V_{base}$ que se calcula como $3.\dfrac{C_{acido}.V_{acido}}{C_{base}}$

A partir de eso se generan los siguientes vectores:

$\begin{align}

	\mathbf{V_{ag, base}} \rightarrow &  \ \{x \in  \mathbb{R} | 0 \leq x \leq V_{base} \}
	
 \\
 \\

	\mathbf{C_{base}} \rightarrow & \ \dfrac{C_{base} * \mathbf{V_{ag, base}} }{\mathbf{V_{ag, base}} + V_{acido} } 

\\

\\


	\mathbf{C_{acido}} \rightarrow & \ \dfrac{C_{acido} * V_{acido} }{\mathbf{V_{ag, base}} + V_{acido} } 

\end{align}$

El pH se calcula a partir de vectorizar el uso de la función `pHfast` de la lib sobre $\mathbf{C_{base}}$ y $\mathbf{C_{acido}}$, en este caso representadas por las variables `conc_na_erlen` y `conc_ac_erlen`. El código redefine al NaOH con cada una de las concentraciónes contenidas en $\mathbf{C_{base}}$.


```julia
	pH=map(conc_na_erlen, conc_ac_erlen) do Cₙ,Cₜ
		especie.conc=Cₜ
		System(
			Neutral(1, Cₙ), #Se define un 'objeto' Neutral, representando al NaOH
			especie # La concentración del Ácido no cambia
			especie
		) |> s -> pHfast(s,.1)
	end
```

\

### Punto de equivalencia

Para calcular el/los volumenes de punto de equivalencia se parte de:

$$C_{acido}.V_{acido} = C_{base}.V_{base}$$
$$\dfrac{C_{acido}.V_{acido}}{C_{base}} = V_{base}$$

Luego, como podemos tener mas de un punto de equivalencia en caso de tener ácidos polipróticos, y sabemos que el volumen del segundo punto de equivalencia deberia corresponder al doble del primero y asi en adelante, definimos el siguiente vector:


$\begin{align}

	\mathbf{Chr} \rightarrow &  \ \{x \in  \mathbb{Z} | 1 \leq x \leq \ce{nH+} \},  \textbf{Chr}\text{ por }\textit{Charge}
	
 \\


 
 \\

\mathbf{V_{base}} =& \dfrac{C_{acido}.V_{acido}.\mathbf{Chr} }{C_{base}}

\end{align}$

Ahora entonces el volumen es un vector, no un escalar.



El pH se calcula de forma similar al caso anterior

```julia

		pH = map(1:times) do t  # t representaría cada uno de los componentes de *Chr*
			println(t)
			acido.conc = (conc_ac * vol_ac )/ (vol[t] + vol_ac) #a)
			base=Neutral(1, t* ((conc_na*vol[t])/(vol[t] + vol_ac))  ) #b)
			System(acido, base) |> s -> pHfast(s,.1)
		end
		(
			vol=vol,
		 	pH=pH
		)
	end
end;
#en a) y b) se tiene en cuenta la dilución.
```




"""


# ╔═╡ f56b3295-d66c-4fbd-ba83-4c54e9d6e871
md"""


Sustancia a titular $(@bind esp Select(especies|> keys |> collect))

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
	xlabel!("Volumen NaOH")
	ylabel!("pH")
end

# ╔═╡ Cell order:
# ╟─4969064a-6f6d-11ec-1d9e-db3f5da90e7e
# ╟─4de35c9b-29bb-4e62-88d0-3e426df2aea5
# ╟─971751a5-b338-4e2a-80bd-c59e76ff3734
# ╟─ca743222-5da5-46d3-bfcd-53fab31acd8c
# ╟─5052e33a-687b-43e4-8437-b8fe392cea87
# ╟─f56b3295-d66c-4fbd-ba83-4c54e9d6e871
# ╟─36d23bf8-99ce-41f9-84d0-d59da6ab44b5
# ╟─e9ae7938-e4cf-4a28-932e-5b956436b53d
