module ParseData
export parse_types, get_sust
	macro parse_types( dict ) 
			dict2 = deepcopy(eval(dict))
			typ = eval(pop!(dict2, :cat))
			esc(typ(values(dict2)...)
	end


	function get_sust(sust, data)
		data[sust]
	end
end





