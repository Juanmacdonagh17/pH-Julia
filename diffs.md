diff --git a/notebooks/pHfunhouse.jl b/notebooks/pHfunhouse.jl
index bf88e75..6cffb98 100644
--- a/notebooks/pHfunhouse.jl
+++ b/notebooks/pHfunhouse.jl
@@ -106,12 +106,6 @@ Concentración molar del Ácido: $(@bind y Slider(0.1: 0.05 :1, show_value = tru
 Concentración molar de la Base: $(@bind z Slider(0.1: 0.05 :1, show_value = true))\
 """
 
-# ╔═╡ 6f8cf5e1-ef18-49b1-93f9-ea2c5b97879c
-md"""
-Mostrar punto de equivalencia: (nada yet)\
-$(@bind eq Select(["Si", "No"]))
-"""
-
 # ╔═╡ e803cbfc-b47c-4168-9ca7-19656e7f7202
 begin
 	bureta = collect(1:0.2:w)
@@ -140,39 +134,19 @@ begin
 			append!(pH_list,(14+log10(-con_list[i])))
 		end
 	end
-	#Punto fantasma en el medio de la titulación:
-	#for i in range(1,length=length(pH_list))
-		#if pH_list[i] == -Inf
-			#pH_list_lie = replace(pH_list,pH_list[i] =>(pH_list[i-1]+pH_list[i+1])/2)
-
-		#end
-	#end
+	for i in range(1,length=length(pH_list))
+		if pH_list[i] == -Inf
+			pH_list_lie = replace(pH_list,pH_list[i] =>(pH_list[i-1]+pH_list[i+1])/2)
+			#pH_list[i] == (pH_list[i-1]+pH_list[i+1])/2
+		end
+	end
 end
 
 # ╔═╡ ab6a95fe-05ff-48f4-a75f-7ae31e2862c5
 scatter(bureta,pH_list, title = "Curva de pH", label = ["pH" "pH"], xlabel = "ml NaOH", ylabel = "pH", legend = false, mode="lines");
 
 # ╔═╡ 445a3f83-f82b-4175-a85c-f84c3bb09d3a
-begin 
-	Plots.plot!((bureta,pH_list), title= "Curva de pH", legend =false);
-	if last(pH_list) > 7
-		if -Inf in pH_list 
-			vol_eq = bureta[findall(isequal(-Inf), pH_list)]
-			Plots.plot!((bureta,pH_list), title= "Curva de pH", legend =false);
-			Plots.plot!((vol_eq,7), title= "Curva de pH", legend =false, markercolor = :red, markershape = :pentagon)
-			annotate!((vol_eq.+4), 7.2, "pH = 7", :color);
-		else
-			v1 = bureta[findall(isequal(last((pH_list[pH_list .< 3]))),pH_list)]
-			v2 = bureta[findall(isequal(first((pH_list[pH_list .> 3]))),pH_list)]
-			vol_eq = (v1.+v2)./2
-			Plots.plot!((bureta,pH_list), title= "Curva de pH", legend =false)
-			Plots.plot!((vol_eq,7), title= "Curva de pH", legend =false, markercolor = :red, markershape = :pentagon)
-			annotate!((vol_eq.+4), 7.2, "pH = 7", :color);
-		end
-	else
-		Plots.plot!((bureta,pH_list), title= "Curva de pH", legend =false)
-	end
-end
+Plots.plot!((bureta,pH_list), title= "Curva de pH", legend =false)
 
 # ╔═╡ 3d7010d3-3ced-4fe6-bcce-01e41d656c57
 #if pH_list[46] == -Inf
@@ -203,7 +177,7 @@ Dict(
 		:viraje=> [10,8.2],
 		:colores=>[:pink,:whitesmoke]
 ),
-	"Naranja de Metilo"=>
+	"Naranja de metilo"=>
 	Dict(
 		:viraje=>[4.4,3.1],
 		:colores=>[:red,:yellow]
@@ -213,12 +187,12 @@ Dict(
 		:viraje=>[12.4,11],
 		:colores=>[:red,:yellow]
 	),
-	"Rojo de Cresol"=>
+	"Rojo de cresol"=>
 	Dict(
 		:viraje=>[8.8,7],
 		:colores=>[:yellow,:red]
 	),
-	"Violeta de Metilo"=>
+	"Violeta de metilo"=>
 	Dict(
 		:viraje=>[1.6,.2],
 		:colores=>[:blueviolet,:yellow]
@@ -1139,14 +1113,8 @@ version = "0.9.1+5"
 # ╟─af14eab0-655e-425b-a483-dc63fc501045
 # ╟─8b512d53-1a96-40f8-88d6-9caf28bbe9ef
 # ╟─1aefcb86-c78f-4bee-9205-75c2b1d972d3
-<<<<<<< HEAD
-# ╟─6f8cf5e1-ef18-49b1-93f9-ea2c5b97879c
-# ╟─e803cbfc-b47c-4168-9ca7-19656e7f7202
-# ╟─8f64a5cd-c290-4972-8714-f914d730f700
-=======
 # ╠═e803cbfc-b47c-4168-9ca7-19656e7f7202
 # ╠═8f64a5cd-c290-4972-8714-f914d730f700
->>>>>>> f608a98287a1017c593920e53ee0d44f75998cce
 # ╟─ab6a95fe-05ff-48f4-a75f-7ae31e2862c5
 # ╟─445a3f83-f82b-4175-a85c-f84c3bb09d3a
 # ╟─3d7010d3-3ced-4fe6-bcce-01e41d656c57
@@ -1154,7 +1122,7 @@ version = "0.9.1+5"
 # ╠═537ea4b1-8cd1-4897-97be-a859d6519b68
 # ╠═43bc566e-13df-4bdf-9933-ae8abf2eda56
 # ╟─1cc940f7-c4a9-4ccc-8c40-94b6694e7fa0
-# ╟─c05192ea-46d5-4a4d-8d7b-aabd5035c241
+# ╠═c05192ea-46d5-4a4d-8d7b-aabd5035c241
 # ╟─e7a769e6-7cbd-4aa3-af4f-eea1dbb34d45
 # ╠═c613f18b-ce91-478d-b8f7-858934f5f932
 # ╟─00000000-0000-0000-0000-000000000001
diff --git a/notebooks/pHfunhouse.jl b/notebooks/pHfunhouse.jl
index bf88e75..6cffb98 100644
--- a/notebooks/pHfunhouse.jl
+++ b/notebooks/pHfunhouse.jl
@@ -106,12 +106,6 @@ Concentración molar del Ácido: $(@bind y Slider(0.1: 0.05 :1, show_value = tru
 Concentración molar de la Base: $(@bind z Slider(0.1: 0.05 :1, show_value = true))\
 """
 
-# ╔═╡ 6f8cf5e1-ef18-49b1-93f9-ea2c5b97879c
-md"""
-Mostrar punto de equivalencia: (nada yet)\
-$(@bind eq Select(["Si", "No"]))
-"""
-
 # ╔═╡ e803cbfc-b47c-4168-9ca7-19656e7f7202
 begin
 	bureta = collect(1:0.2:w)
@@ -140,39 +134,19 @@ begin
 			append!(pH_list,(14+log10(-con_list[i])))
 		end
 	end
-	#Punto fantasma en el medio de la titulación:
-	#for i in range(1,length=length(pH_list))
-		#if pH_list[i] == -Inf
-			#pH_list_lie = replace(pH_list,pH_list[i] =>(pH_list[i-1]+pH_list[i+1])/2)
-
-		#end
-	#end
+	for i in range(1,length=length(pH_list))
+		if pH_list[i] == -Inf
+			pH_list_lie = replace(pH_list,pH_list[i] =>(pH_list[i-1]+pH_list[i+1])/2)
+			#pH_list[i] == (pH_list[i-1]+pH_list[i+1])/2
+		end
+	end
 end
 
 # ╔═╡ ab6a95fe-05ff-48f4-a75f-7ae31e2862c5
 scatter(bureta,pH_list, title = "Curva de pH", label = ["pH" "pH"], xlabel = "ml NaOH", ylabel = "pH", legend = false, mode="lines");
 
 # ╔═╡ 445a3f83-f82b-4175-a85c-f84c3bb09d3a
-begin 
-	Plots.plot!((bureta,pH_list), title= "Curva de pH", legend =false);
-	if last(pH_list) > 7
-		if -Inf in pH_list 
-			vol_eq = bureta[findall(isequal(-Inf), pH_list)]
-			Plots.plot!((bureta,pH_list), title= "Curva de pH", legend =false);
-			Plots.plot!((vol_eq,7), title= "Curva de pH", legend =false, markercolor = :red, markershape = :pentagon)
-			annotate!((vol_eq.+4), 7.2, "pH = 7", :color);
-		else
-			v1 = bureta[findall(isequal(last((pH_list[pH_list .< 3]))),pH_list)]
-			v2 = bureta[findall(isequal(first((pH_list[pH_list .> 3]))),pH_list)]
-			vol_eq = (v1.+v2)./2
-			Plots.plot!((bureta,pH_list), title= "Curva de pH", legend =false)
-			Plots.plot!((vol_eq,7), title= "Curva de pH", legend =false, markercolor = :red, markershape = :pentagon)
-			annotate!((vol_eq.+4), 7.2, "pH = 7", :color);
-		end
-	else
-		Plots.plot!((bureta,pH_list), title= "Curva de pH", legend =false)
-	end
-end
+Plots.plot!((bureta,pH_list), title= "Curva de pH", legend =false)
 
 # ╔═╡ 3d7010d3-3ced-4fe6-bcce-01e41d656c57
 #if pH_list[46] == -Inf
@@ -203,7 +177,7 @@ Dict(
 		:viraje=> [10,8.2],
 		:colores=>[:pink,:whitesmoke]
 ),
-	"Naranja de Metilo"=>
+	"Naranja de metilo"=>
 	Dict(
 		:viraje=>[4.4,3.1],
 		:colores=>[:red,:yellow]
@@ -213,12 +187,12 @@ Dict(
 		:viraje=>[12.4,11],
 		:colores=>[:red,:yellow]
 	),
-	"Rojo de Cresol"=>
+	"Rojo de cresol"=>
 	Dict(
 		:viraje=>[8.8,7],
 		:colores=>[:yellow,:red]
 	),
-	"Violeta de Metilo"=>
+	"Violeta de metilo"=>
 	Dict(
 		:viraje=>[1.6,.2],
 		:colores=>[:blueviolet,:yellow]
@@ -1139,14 +1113,8 @@ version = "0.9.1+5"
 # ╟─af14eab0-655e-425b-a483-dc63fc501045
 # ╟─8b512d53-1a96-40f8-88d6-9caf28bbe9ef
 # ╟─1aefcb86-c78f-4bee-9205-75c2b1d972d3
-<<<<<<< HEAD
-# ╟─6f8cf5e1-ef18-49b1-93f9-ea2c5b97879c
-# ╟─e803cbfc-b47c-4168-9ca7-19656e7f7202
-# ╟─8f64a5cd-c290-4972-8714-f914d730f700
-=======
 # ╠═e803cbfc-b47c-4168-9ca7-19656e7f7202
 # ╠═8f64a5cd-c290-4972-8714-f914d730f700
->>>>>>> f608a98287a1017c593920e53ee0d44f75998cce
 # ╟─ab6a95fe-05ff-48f4-a75f-7ae31e2862c5
 # ╟─445a3f83-f82b-4175-a85c-f84c3bb09d3a
 # ╟─3d7010d3-3ced-4fe6-bcce-01e41d656c57
@@ -1154,7 +1122,7 @@ version = "0.9.1+5"
 # ╠═537ea4b1-8cd1-4897-97be-a859d6519b68
 # ╠═43bc566e-13df-4bdf-9933-ae8abf2eda56
 # ╟─1cc940f7-c4a9-4ccc-8c40-94b6694e7fa0
-# ╟─c05192ea-46d5-4a4d-8d7b-aabd5035c241
+# ╠═c05192ea-46d5-4a4d-8d7b-aabd5035c241
 # ╟─e7a769e6-7cbd-4aa3-af4f-eea1dbb34d45
 # ╠═c613f18b-ce91-478d-b8f7-858934f5f932
 # ╟─00000000-0000-0000-0000-000000000001
