### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 99cb7cf8-a9df-11ef-221d-87167b2e404a
begin
	using PlutoUI,DataFrames, CSV, Random, PrettyTables, Format, JuMP, HiGHS, HypertextLiteral
	Print()
end

# ╔═╡ d5be6e8a-6e84-4a33-87a9-7b972bf92edc
@bind go Button("Recompute")

# ╔═╡ 03bffb9e-d5a2-458e-ad2b-cbc87002bc82
begin
	function makeSampleRotationOrderBook()
	    CSV.read(download("https://raw.githubusercontent.com/adannenberg/PlutoExport/main/spx_constituents.csv"), DataFrame)
	    nsec = length(tickers.ticker)
	    npm = 9
	    numOrders = npm * nsec
	    #Random.seed!(456) #321
	    v = Float64.(rand([200, 700, 1200], numOrders))
	    p = rand(1:nsec, numOrders)
	    orderBook = DataFrame(
	        pm=string.(rand(1:npm, numOrders)),
	        order_id=string.(1:numOrders),
	        sell_ticker=tickers.ticker[p],
	        sell_mv=v,
	        buy_ticker=tickers.ticker[vcat(p[end], p[1:end-1])],
	        buy_mv=v
	    )
	    orderBook.sell_mv[orderBook.pm.=="1"] .*= 2
	    orderBook.buy_mv[orderBook.pm.=="2"] .*= 3
	    # remove rows in which buy_ticker = sell_ticker
	    indx = findall(orderBook.buy_ticker .== orderBook.sell_ticker)
	    if length(indx) > 0
	        orderBook = orderBook[Not(indx), :]
	    end
	
	    return orderBook
	end
	
	function convertRotationOrderBookToRegularOrderBook(rotation_ob::DataFrame)
	    regular_ob = DataFrame(pm=String[], order_id=String[], side=String[], ticker=String[], mv=Float64[])
	    for i in 1:nrow(rotation_ob)
	        push!(regular_ob, (rotation_ob.pm[i], "buy_" * string(rotation_ob.order_id[i]), "buy", rotation_ob.buy_ticker[i], rotation_ob.buy_mv[i]))
	        push!(regular_ob, (rotation_ob.pm[i], "sell_" * string(rotation_ob.order_id[i]), "sell", rotation_ob.sell_ticker[i], rotation_ob.sell_mv[i]))
	    end
	    return regular_ob
	end
	
	function ob_flattener(orderBook::DataFrame)
	    ob = copy(orderBook)
	    
	    if "order_id" ∈ names(ob)
	        ob = ob[:,Not(:order_id)]
	    end
	    ob.sign .= ifelse.(ob.side .== "buy", 1, -1)
	    pms = sort(unique(ob.pm))
	    tickers = sort(unique(ob.ticker))
	    ob2 = copy(ob[1:-1,:]) #just an initialization trick
	    ob2 = ob2[:,Not(:sign)]
	    for pm in pms, ticker in tickers
	       indx = findall(ob.pm .== pm .&& ob.ticker .== ticker)
	       if length(indx) > 0
	            mv = sum(ob.mv[indx] .* ob.sign[indx])
	            if mv > 1
	                push!(ob2, (pm = pm, side = "buy", ticker = ticker, mv = mv))
	            elseif mv < -1
	                push!(ob2, (pm = pm, side = "sell", ticker = ticker, mv = -mv))
	            end
	        end
	    end
	    return ob2
	end
	
	function punchline(ob::DataFrame)
	    pms = sort(unique(ob.pm))
	    pm_info = DataFrame(pm=pms, order_gmv = zeros(length(pms)), executed_gmv=zeros(length(pms)), order_nmv = zeros(length(pms)), executed_nmv=zeros(length(pms)))
	    for i in 1:nrow(pm_info)
	        pm = pm_info.pm[i]
	        pm_indx = findall(ob.pm .== pm)
	        if length(pm_indx) > 0
	            pm_info.order_gmv[i] = sum(ob.mv[pm_indx])
	            pm_info.executed_gmv[i] = sum(abs.(ob.execmv[pm_indx]))
	            pm_info.order_nmv[i] = sum(ob.signedmv[pm_indx])
	            pm_info.executed_nmv[i] = sum(ob.execmv[pm_indx])
	        end
	    end
	    push!(pm_info, (pm="total", order_gmv = sum(pm_info.order_gmv), executed_gmv=sum(pm_info.executed_gmv), order_nmv = sum(pm_info.order_nmv), executed_nmv=sum(pm_info.executed_nmv)))
	    #pm_info.xnmv2xgmv = string.(round.(100 * pm_info.executed_nmv ./ pm_info.executed_gmv, digits=1)) .* "%"
	    #pm_info.xgmv2ogmv = string.(round.(100*pm_info.executed_gmv ./ pm_info.order_gmv, digits=1)) .* "%"
		pm_info[:,"executed nmv / gmv"] = string.(round.(100 * pm_info.executed_nmv ./ pm_info.executed_gmv, digits=1)) .* "%"
	    pm_info[:,"executed gmv / order gmv"] = string.(round.(100*pm_info.executed_gmv ./ pm_info.order_gmv, digits=1)) .* "%"
		rename!(pm_info, ["pm","order gmv", "executed gmv", "order nmv", "executed nmv", "executed nmv / gmv", "executed gmv / order gmv"])

	    for c in 2:5
	        pm_info[:,c] = Int64.(round.(pm_info[:,c]))
	    end
	    return pm_info
	end
	
	function tMat2ExecutionReport(tMat::Matrix, pms::Vector{String}, tickers::Vector{String})
	    #to be called within crossing_engine so ob, tickers, pms will be within scope
	    # creates an execution report of internalized trades with format DataFrame(mv=Float64[], ticker=String[], buying_pm=String[], selling_pm=String[])
	    n_ticker = length(tickers)
	    executions = DataFrame(mv=Float64[], ticker=String[], buying_pm=String[], selling_pm=String[])
	    for c in 1:n_ticker
	        ticker = tickers[c]
	        nmvs = copy(tMat[:,c])
	        gmvs = abs.(nmvs)
	        while true
	            buy_indx = findall(nmvs .> 0)
	            sell_indx = findall(nmvs .< 0)
	            if length(buy_indx) * length(sell_indx) == 0
	                break
	            end
	            mv = min(nmvs[buy_indx[1]], -nmvs[sell_indx[1]])
	            buying_pm = pms[buy_indx[1]]
	            selling_pm = pms[sell_indx[1]]
	            push!(executions, (mv = mv, ticker = ticker, buying_pm = buying_pm, selling_pm = selling_pm))
	            nmvs[buy_indx[1]] -= mv
	            nmvs[sell_indx[1]] += mv
	            gmvs = abs.(nmvs)
	        end
	    end
	    return executions
	end
	
	function crossing_engine(ob::DataFrame)
	    x = copy(ob)
	    tickers = sort(unique(x.ticker))
	    pms = sort(unique(x.pm))
	    n_pm = length(pms)
	    n_ticker = length(tickers)
	    oMat = zeros(n_pm, n_ticker)
	    for i in 1:nrow(x)
	        r = findall(pms .== x.pm[i])[1]
	        c = findall(tickers .== x.ticker[i])[1]
	        x.side[i] == "buy" ? sgn = 1 : sgn = -1
	        oMat[r, c] = sgn * x.mv[i]
	    end
	
	    sMat = sign.(oMat)
	
	    model = Model(HiGHS.Optimizer)
	
	    @variable(model, tmv[1:n_pm, 1:n_ticker])
	    
	    # order-level restrictions
	    zindx = findall(oMat .== 0)
	    if length(zindx) > 0
	        @constraint(model, c0, tmv[zindx] .== 0)
	    end
	    nzindx = findall(abs.(oMat) .> 0)
	    if length(nzindx) > 0
	        @constraint(model, c3, tmv[nzindx] .* sMat[nzindx] .>= 0)  #ensures that all trades are of the right sign
	        @constraint(model, c4, tmv[nzindx] .* sMat[nzindx] .<= oMat[nzindx] .* sMat[nzindx]) # ensures that trades are not bigger than orders
	    end
	    # pm-level restrictions
	    indx = findall(sum(oMat, dims=2)[:] .== 0)
	    if length(indx) > 0
	        # c1a makes each each pm with a mn list have mn execns
	        @constraint(model, c1a, sum(tmv, dims=2)[indx, 1] .== 0)
	    end
	    
	    indx = findall(sum(oMat, dims=2)[:] .!= 0)
	    if length(indx) > 0
	        # c1b makes each pm with an outright list have nmv ∈ [-Inf, +nmv_target] for nmv_target > 0, nmv ∈ [nmv_target, +Inf] for nmv_target < 0
	        @constraint(model, c1b, sum(tmv, dims=2)[indx,1] ./ sum(oMat, dims=2)[indx,1] .<= 1) 
	        # c1c complements c1b and makes each pm have sign (nmv) = sign(nmv_target)
	        # in conjunction with c1b we now have the requirement that nmv lie between 0 and nmv_target
	        @constraint(model, c1c, sum(tmv, dims=2)[indx,1] .* sum(oMat, dims=2)[indx,1] .>= 0)
	    end
	    # by requiring that total trades in each ticker = 0 we ensure all trading is internalized
	    @constraint(model, c2, sum(tmv, dims=1) .== 0)
	    # objective is to maximize gmv traded subject to constraints
	    @objective(model, Max, sum(sMat .* tmv))
	    set_silent(model)
	    optimize!(model)
	    tMat = value.(tmv) #matrix of trade match_values
	    # now let's take mn_ob and add a net executed value column from tMat
	
	    x.sign .= ifelse.(x.side .== "buy", 1, -1)
	    x.signedmv = x.mv .* x.sign
	
	    x.execmv = zeros(nrow(x))
	    for i in 1:n_pm, j in 1:n_ticker
	        pm = pms[i]
	        ticker = tickers[j]
	        if tMat[i,j] != 0
	            r = findall(x.pm .== pm .&& x.ticker .== ticker)
	            if length(r) != 1
	                println("PROBLEM: for pm = $pm, ticker = $ticker optimzn gives tradevalue $(tMat[i,j]) but there's no corresponding order")
	            else
	                x.execmv[r[1]] = tMat[i,j]
	            end
	        end
	    end
	    x.opensignedmv = x.signedmv .- x.execmv
	    executionReport = tMat2ExecutionReport(tMat, pms, tickers)
	    return x, tMat, executionReport
	end

	function ptf(v,i,j)
	    if !(v isa String)
	        return replace(string(Int(round(v))), r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
	    else
	        return v
	    end
	end

	function display_df_scrollable(df::DataFrame; 
	                             title::String = "my title",
	                             max_height::Int=400,
	                             max_rows::Union{Int,Nothing}=nothing,
								 showRowLabels::Bool = true,
	                             scrollable::Bool=true)
	    
	    # If max_rows is specified, limit the displayed rows
	    display_df = isnothing(max_rows) ? df : first(df, max_rows)
	    
	    # Get the basic table HTML
		if showRowLabels == true
		    table_html = pretty_table(String, display_df, 
		        backend=Val(:html),
		        formatters=ptf,
		        title = title,
		        title_alignment = :C,
		        row_labels = 1:nrow(display_df),
		        show_subheader = false,
		        standalone=false,
		        tf=tf_html_default)
		else
			table_html = pretty_table(String, display_df, 
		        backend=Val(:html),
		        formatters=ptf,
		        title = title,
		        title_alignment = :C,
		        show_subheader = false,
		        standalone=false,
		        tf=tf_html_default)
		end
	    
	    # Replace the existing table tags with ones that have our styles, removing borders
	    styled_table = replace(table_html, 
	        "<table>" => """<table style="width:100%; border-collapse:collapse; font-size:16px; border-spacing:0; position:relative;  ">""") #border:none;
	    
	    # Adjust header styling based on scrollable parameter, removing borders
	    header_top = scrollable ? "top:40px;" : "top:0;"
	    styled_table = replace(styled_table,
	        "<th" => """<th style="position:sticky; $(header_top) background:white; padding:8px; text-align:center; font-size:16px; " """) #border:none;
	    
	    styled_table = replace(styled_table,
	        "<td" => """<td style="padding:8px; text-align:center; font-size:18px; " """) #border:none;
	    
	    # Adjust caption styling based on scrollable parameter
	    caption_position = scrollable ? "position: sticky; top: 0;" : ""
	    styled_table = replace(styled_table,
	        """<caption style = "text-align: center;">""" => 
	        """<caption style="text-align: center; font-family: Georgia, serif; font-size: 24px; font-weight: bold; color: #2c3e50; padding: 10px; $(caption_position) background: white; z-index: 1; ">""") #border:none;
	    
	    # Wrap in appropriate div without borders
	    if scrollable
	        return @htl("""
	            <div style="max-height:$(max_height)px; overflow-y:auto; display:block;">
	                $(HTML(styled_table))
	            </div>
	        """)
	    else
	        return @htl("""
	            <div style="display:block;">
	                $(HTML(styled_table))
	            </div>
	        """)
	    end
	end
		
	Print()
end

# ╔═╡ 77adcf13-1ba4-4e6a-8744-3636d0f67ca1
begin
	go
	ob = convertRotationOrderBookToRegularOrderBook(makeSampleRotationOrderBook())
	ob = ob_flattener(ob)
	ob_exec, tMat, executionReport = crossing_engine(ob);
	Print()
end

# ╔═╡ 61e99319-3be7-4646-82bc-db7e34f50bdc
display_df_scrollable(punchline(ob_exec),title="summary statistics of internal cross",scrollable=false,showRowLabels=false)

# ╔═╡ 60422a72-068a-46e2-84d8-f6900a67d73e
display_df_scrollable(ob, title = "simulated order book of $(format(nrow(ob),commas=true)) orders among $(length(unique(ob.pm))) PMs", showRowLabels = true)

# ╔═╡ 6003dec6-7270-47e1-833d-3f283a063181
display_df_scrollable(executionReport, title = "internal crosses", showRowLabels = true)

# ╔═╡ 07e3477a-c0df-499a-99f3-bb2eff71ec5b
begin
	indx = findall(ob_exec.opensignedmv .!= 0)
	x = ob_exec[indx, [:pm, :side, :ticker, :opensignedmv]]
	x.opensignedmv = abs.(x.opensignedmv)
	rename!(x, [:pm,:side,:ticker,:mv])
	display_df_scrollable(x, title = "remaining open orders")
end

# ╔═╡ 492db5f9-6916-41e4-9088-d9b31de42996
html"""<style>
	main {
    	max-width: 1000px !important;
    	margin-left: 100px !important;
	}
	"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Format = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
HiGHS = "87dc4568-4c63-4d18-b0c0-bb2238e4078b"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
CSV = "~0.10.15"
DataFrames = "~1.7.0"
Format = "~1.3.7"
HiGHS = "~1.12.1"
HypertextLiteral = "~0.9.5"
JuMP = "~1.23.5"
PlutoUI = "~0.7.60"
PrettyTables = "~2.4.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.2"
manifest_format = "2.0"
project_hash = "7187af67ecd89b50b5516617312fa1d678416b7c"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "f1dff6729bc61f4d49e140da1af55dcd1ac97b2f"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.5.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8873e196c2eb87962a2048b3b8e08946535864a1"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+2"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "TranscodingStreams"]
git-tree-sha1 = "e7c529cc31bb85b97631b922fa2e6baf246f5905"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.8.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "7878ff7172a8e6beedd1dea14bd27c3c6340d361"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.22"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

    [deps.ForwardDiff.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.HiGHS]]
deps = ["HiGHS_jll", "MathOptInterface", "PrecompileTools", "SparseArrays"]
git-tree-sha1 = "02fad6652a24cd1356c5dc000c3cca297e13482a"
uuid = "87dc4568-4c63-4d18-b0c0-bb2238e4078b"
version = "1.12.1"

[[deps.HiGHS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "cc963ae42b15ccd2536fb7a6254c0328b404347c"
uuid = "8fd58aa0-07eb-5a78-9b36-339c94fd15ea"
version = "1.8.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InlineStrings]]
git-tree-sha1 = "45521d31238e87ee9f9732561bfee12d4eebd52d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.2"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuMP]]
deps = ["LinearAlgebra", "MacroTools", "MathOptInterface", "MutableArithmetics", "OrderedCollections", "PrecompileTools", "Printf", "SparseArrays"]
git-tree-sha1 = "866dd0bf0474f0d5527c2765c71889762ba90a27"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "1.23.5"

    [deps.JuMP.extensions]
    JuMPDimensionalDataExt = "DimensionalData"

    [deps.JuMP.weakdeps]
    DimensionalData = "0703355e-b756-11e9-17c0-8b28908087d0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "PrecompileTools", "Printf", "SparseArrays", "SpecialFunctions", "Test", "Unicode"]
git-tree-sha1 = "e065ca5234f53fd6f920efaee4940627ad991fb4"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.34.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "a2710df6b0931f987530f59427441b21245d8f5e"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Profile]]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "d0553ce4031a081cc42387a9b9c8441b7d99f32d"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.7"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a6b1675a536c5ad1a60e5a5153e1fee12eb146e3"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─d5be6e8a-6e84-4a33-87a9-7b972bf92edc
# ╟─61e99319-3be7-4646-82bc-db7e34f50bdc
# ╟─60422a72-068a-46e2-84d8-f6900a67d73e
# ╟─6003dec6-7270-47e1-833d-3f283a063181
# ╟─07e3477a-c0df-499a-99f3-bb2eff71ec5b
# ╟─77adcf13-1ba4-4e6a-8744-3636d0f67ca1
# ╟─03bffb9e-d5a2-458e-ad2b-cbc87002bc82
# ╟─492db5f9-6916-41e4-9088-d9b31de42996
# ╟─99cb7cf8-a9df-11ef-221d-87167b2e404a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
