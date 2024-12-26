### A Pluto.jl notebook ###
# v0.20.4

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

# ╔═╡ 86d497c4-109a-4642-996d-50fb9bb6c6d0
begin
	  if Sys.isapple()
		push!(LOAD_PATH, "/Users/alex/Library/CloudStorage/Dropbox/Code/Julia/Projects/")
     elseif Sys.islinux()
  		 push!(LOAD_PATH,"/home/alex/Dropbox/Code/Julia/Projects/")
     end
	import Pkg
    Pkg.activate()
	using Plots, PlutoUI, BenchmarkTools, DataFrames, Dates,  Revise, ahd_funcs
	Print()
end

# ╔═╡ a784bc70-4b47-4154-a2c2-56b676ed91e8
html"""
<style>
div.edit_or_run {
    display: none !important;
    visibility: hidden !important;
    width: 0 !important;
    height: 0 !important;
    margin: 0 !important;
    padding: 0 !important;
    border: none !important;
    pointer-events: none !important;
    position: absolute !important;
    clip: rect(0, 0, 0, 0) !important;
    overflow: hidden !important;
}
</style>
"""

# ╔═╡ 66f07370-1190-4f2a-a509-fb0bdc72653d
md"""
## Conversion Calculator w/ Rate Inputs
"""

# ╔═╡ 707e1bdc-e84b-4d5c-b44a-a6a61b877e85
PlutoUI.combine() do Child
md"""
**Implied Dividend Section**:

Compute implied dividend(s) that reprice market conversion price:$(@bind doImpDivs CheckBox(default = false))
"""
end

# ╔═╡ 877df67b-6950-4a5c-881c-b252c77e2707
if doImpDivs
	@bind maxDiv2Px confirm(
		PlutoUI.combine() do Child
			md"""
			Max value of divamt / px to consider:
			$(Child(Slider(0:.01:.9, show_value = true, default = .1)))
			"""
		end
	)
	#maxDiv2Px = maxDiv2Px[1]
end

# ╔═╡ c1003ee9-c023-41ea-a4ef-b1493287afa3
begin
	nbsp = html"&nbsp"
	Print()
end

# ╔═╡ bd181597-a3dc-4a15-91f0-9611a330c8db
@bind inputs confirm(
		PlutoUI.combine() do Child
		md"""
		
		**Conversion Parameters**:
		
		Stock Price:  $(Child(NumberField(0:.01:10000, default=100.00)))
		$nbsp Strike Price: $(Child(NumberField(5.0:1000.0,default = 120.00)))
		$nbsp Volatility: $(Child(NumberField(0.0:500.0,default = 25.00)))
		
		Expiration Date: $(Child(Select(bizDayAdjustedThirdFriday.(today()+Month.(0:60)), default = bizDayAdjustedThirdFriday(today() .+ Year(1)))))
		$nbsp 
		Trade date: $(Child(Select(unique(bd.(today() - Day(180):Day(1):today() + Day(180))), default = bd(today()))))
		
		Funding Rate: $(Child(NumberField(0.0:.1:10.0,default = 5.0)))
		$nbsp Repo Rate: $(Child(NumberField(-500.0:.1:10.0, default = -1.0)))
		
		Dividend amounts:$(Child(NumberField(0.0:100.0, default = 2.00)))
		$(Child(NumberField(0.0:100.0, default = 2.00)))
		$(Child(NumberField(0.0:100.0, default = 0.00)))
		
		Ex-dividend Dates: 
		$(Child(Select(unique(bd.(today()-Day(180):Day(1):today()+Year(5))), default = bd(today()+Month(3)))))
		$nbsp 
		$(Child(Select(unique(bd.(today():Day(1):today()+Year(5))), default = bd(today()+Month(6)))))
		$nbsp 
		$(Child(Select(unique(bd.(today():Day(1):today()+Year(5))), default = bd(today()+Month(9)))))
		
		Dividend Pay Dates:
		$(Child(Select(unique(bd.(today():Day(1):today()+Year(5))), default = bd(today()+Month(3)+Week(1)))))
		$nbsp 
		$(Child(Select(unique(bd.(today():Day(1):today()+Year(5))), default = bd(today()+Month(6)+Week(1)))))
		$nbsp 
		$(Child(Select(unique(bd.(today():Day(1):today()+Year(5))), default = bd(today()+Month(9)+Week(1)))))

		Tree type: $(Child(Select(["CRR","LR","Trinomial"], default = "CRR")))
		$nbsp
		Steps per day: $(Child(NumberField(1.0:50.0,default = 1.0)))

		Market Price of Conversion: $(Child(NumberField(-50.0:.1:50.0, default = 2.0)))
		"""
	end
	)

# ╔═╡ bde7c9cc-a124-419b-ab02-6c2273e87a26
begin
	(S, K, vol_pct, dt_expiry, dt_trade, r_csa_mm_pct, r_repo_mm_pct, divamt1, divamt2, divamt3, dt_exdiv1, dt_exdiv2, dt_exdiv3, dt_divpmt1, dt_divpmt2, dt_divpmt3, treeType, stepsPerDay, mktConversionPx) = inputs

	Print()
end

# ╔═╡ f5c2fd63-5fc6-40bd-9a7f-9731ea5636b8
begin
	divamt = [divamt1; divamt2; divamt3]
	dt_exdiv = [dt_exdiv1; dt_exdiv2; dt_exdiv3]
	dt_divpmt = [dt_divpmt1; dt_divpmt2; dt_divpmt3]
	local idx = findall(divamt .> 0)
	if length(idx) == 0
		divamt = [divamt[1]]
		dt_exdiv = [dt_exdiv[1]]
		dt_divpmt = [dt_divpmt[1]]
	else
		divamt = divamt[idx]
		dt_exdiv = dt_exdiv[idx]
		dt_divpmt = dt_divpmt[idx]
	end
	local indx = findall(dt_exdiv .<= dt_expiry)
	if length(indx) > 0
		divamt = divamt[indx]
		dt_exdiv = dt_exdiv[indx]
		dt_divpmt = dt_divpmt[indx]
	else
		divamt = [0]
		dt_exdiv = [bd(dt_expiry-Day(7))]
		dt_divpmt = [bd(dt_expiry-Day(7))]
	end
	Print()
end

# ╔═╡ d48514f9-8e54-4666-a0ae-4a5a62daae2f
if doImpDivs
	#px = conversion_Basic(S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_exp, stepsPerDay, dt_exdiv, dt_divpmt, divamt, dt_trade, earlyExPut, earlyExCall)[1]
	
	impDiv = implied_dividend_from_conversion_am(mktConversionPx, treeType, S, K, r_csa_mm_pct, r_repo_mm_pct,vol_pct, dt_trade, dt_expiry, dt_exdiv, dt_divpmt, stepsPerDay, maxDiv2Px[1])
	Print("implied dividend = $(round(impDiv,digits=2))")
end

# ╔═╡ 0c482b6d-cd42-4036-8cea-7582c774eef3
if doImpDivs
	divrng = collect(LinRange(0,maxDiv2Px[1],10)) .* S
	#divrng = collect(0:.01:maxDiv2Px[1]) .* S
	divDF = DataFrame(div = divrng, conv = zeros(length(divrng)) .+ NaN)
	ndiv = length(dt_exdiv)
	Threads.@threads for r in 1:nrow(divDF)
		divDF.conv[r] = conversion_am(treeType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divDF.div[r] * ones(ndiv), dt_exdiv, dt_divpmt, stepsPerDay)
	end
	scatter(divDF.div,divDF.conv,xlabel="dividend",ylabel="conv px", label="conv px as function of dividend", color=:blue)
	plot!(divDF.div,divDF.conv,label="",color=:blue)
	hline!([mktConversionPx],color=:red, label="mkt conv px = $mktConversionPx")
end

# ╔═╡ 9f048c1e-8e5c-4e78-a587-3e8566fdc098
md"""
**Misc:**:

Include value of unexercised options in pnl calcs: $(@bind includeUnexOptVals CheckBox(default = true))

Use American pricer for unexercised option values (slower): $(@bind useAmericanPricerForUnexOpts CheckBox(default = false))

Assume early exercise is possible for:

Put holder $(@bind earlyExPut CheckBox(default = true))
$nbsp $nbsp Call holder $(@bind earlyExCall CheckBox(default = true))
"""

# ╔═╡ 3046063a-0396-4e8c-a602-b9f8f17295d4
@bind graph_inputs confirm(
	PlutoUI.combine() do Child
	md"""
	**Plot Selection**:
	
	Plot conversion price vs repo rate:$(Child(CheckBox(default = false)))
	$nbsp
	for repo rng  $(Child(Slider(-1000:0, show_value = true, default = -10)))
	$nbsp
	to  $(Child(Slider(0:100, show_value = true, default = 1)))
	
	Do delta/gamma graphs:$(Child(CheckBox(default = false)))
	
	Do exercise boundary graphs:$(Child(CheckBox(default = false)))

	Do mtm pnl contour graph (not available):$(Child(CheckBox(default = false)))
	
	"""
	end
)

# ╔═╡ 39493c2e-0870-417b-94bc-3658c1c18d34
begin
	(doRepoGraph,minRate,maxRate,doRiskGraphs,doExBoundaryGraphs, doMtmPnlGraph) = graph_inputs
	Print()
end

# ╔═╡ 6413cc25-ac66-4989-8307-312a35243d1f
begin
	returnRangePct = LinRange(minRate,maxRate,2*Threads.nthreads()) # just for delta/gamma graphs
	days = (dt_expiry - dt_trade).value
	r_repo = Rmm2Rcc(r_repo_mm_pct, days) / 100
	r_csa = Rmm2Rcc(r_csa_mm_pct, days) /100
	vol = vol_pct / 100.0

	max_totrtn = -1 + exp((r_repo - 0.5 * vol^2) * (dt_expiry - dt_trade).value / 365 + vol * sqrt((dt_expiry - dt_trade).value / 365) * 3)
	min_totrtn = -1 + exp((r_repo - 0.5 * vol^2) * (dt_expiry - dt_trade).value / 365 + vol * sqrt((dt_expiry - dt_trade).value / 365) * -3)
	max_rtn_pct = 10 * ceil(10 * max_totrtn)
	min_rtn_pct = 10 * floor(10 * min_totrtn)


	numPrices = round(Int, 1 + (max_rtn_pct - min_rtn_pct) / 10) # number of x-axis points for graph
	numDates = maximum([20; round(Int, (dt_expiry - dt_trade).value / 10)]) # number of y-axis points for graph

	#American model results
	conv = conversion_am(treeType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
	
	p = option_am(treeType, "put", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1]
	
	c = option_am(treeType, "call", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1]
	
	function isbad(x)
		if ismissing(x) || isnan(x)
			return true
		else
			return false
		end
	end
	if !isbad(mktConversionPx)
		irr = implied_repo_from_conversion_am(mktConversionPx, treeType, S, K, r_csa_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
	else
		irr = r_repo_mm_pct
		convMkt = conv
	end

	#European model results
	conv_bs = conversion_BS(S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
	
	p_bs = BS("put", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
	
	c_bs = BS("call", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
	
	if !isbad(mktConversionPx)
		irr2 = implied_repo_from_conversion_BS(mktConversionPx, S, K, r_csa_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
	else
		irr2 = r_repo_mm_pct
	end

	#=
	if !ismissing(pExPx)
		pExPxSynth = pExPx
		local indx = findall(pExDt .>= dt_exdiv)
		if length(indx) > 0
			pExPxSynth = pExPxSynth + sum(divamt[indx])
		end
		pExRtn = log(pExPxSynth/S)
		pExTime = (pExDt - dt_trade).value/365
		pExSigs = (pExRtn - (r_repo - vol^2/2) * pExTime) / (vol * sqrt(pExTime))
		if pExRtn > 0
			pExDir = "above"
		else
			pExDir = "below"
		end
		str1 = "earliest put exercise seen by valuation is at " * string(round(pExPx, digits = 2)) * " on " * Dates.format(pExDt,"u d, Y")
		if isfinite(pExSigs)
			str1 = str1 * " which is " * string(abs(round(Int, pExSigs))) * "σ "* pExDir *" spot"
		end
	else
		str1 = " "
	end
	if !ismissing(cExPx)
		cExPxSynth = cExPx
		local indx = findall(cExDt .>= dt_exdiv)
		if length(indx) > 0
			cExPxSynth = cExPxSynth + sum(divamt[indx])
		end
		cExRtn = log(cExPxSynth/S)
		cExTime = (cExDt - dt_trade).value/365
		cExSigs = (cExRtn - (r_repo - vol^2/2) * cExTime) / (vol * sqrt(cExTime))
		if cExRtn > 0
			cExDir = "above"
		else
			cExDir = "below"
		end
		str2 = "earliest call exercise seen by valuation is at " * string(round(cExPx, digits = 2)) * " on " * Dates.format(cExDt,"u d, Y")
		if isfinite(cExSigs)
			str2 = str2 * " which is " * string(abs(round(Int, cExSigs))) * "σ "* cExDir *" spot"
		end
	else
		str2 = " "
	end
	=#
	
	if isbad(mktConversionPx)
		dum = convMkt
	else
		dum = mktConversionPx
	end
	message = string(
		"** American model **\n",
		"With user-input repo rate of ", round(r_repo_mm_pct,digits=2),"%,",
		"\nconversion = ",
		round(conv, digits = 3),
		"\nput = ",
		round(p, digits = 2),
		"\ncall = ",
		round(c, digits = 2),
		"\n\nBut note that user-input conversion price of ", mktConversionPx," implies a repo rate of ", round(irr,digits = 2),"%",
		
	)
	#=,
		"\n" * str1,
		"\nexpected put life is " * string(round(putTau)) * " days, i.e. until " * Dates.format(dt_trade + Day(round(putTau)),"Y-m-d"),
		"\n" * str2,
		"\nexpected call life is " * string(round(callTau)) * " days, i.e. until " * Dates.format(dt_trade + Day(round(callTau)),"Y-m-d"),
	)
	=#
	message = string(
		message,
		"\n\n** European (BS) model **\n",
		"With user-input repo rate of ", round(r_repo_mm_pct,digits=2),"%,",
		"\nconversion = ",
		round(conv_bs, digits = 3),
		"\nput = ",
		round(p_bs, digits = 2),
		"\ncall = ",
		round(c_bs, digits = 2),
		"\n\nBut note that user-input conversion price of ", mktConversionPx," implies a repo rate of ", round(irr2,digits = 2),"%",
		"\n(using BS and, therefore, ignoring early exercise)"
		
	)
	#=
	message = string(
		message,
		"\n",
		"\nusing European model:",
		"\nconversion = ",
		round(conv_bs, digits = 4),
		"\nput = ",
		round(p_bs, digits = 2),
		"\ncall = ",
		round(c_bs, digits = 2),
		"\nimplied repo rate from mkt conversion price is ",
		round(irr2, digits = 2),
		"%\n",
	)
	=#
	Print()
end

# ╔═╡ 03142ff5-32f4-44a7-a11f-38d55ce06712
Print(message)

# ╔═╡ 97498b0f-b148-41e7-b16e-abb6447f0b61
if doRepoGraph
	r_repo_pct_rng = LinRange(minRate,maxRate,30)	
	conversions = zeros(length(r_repo_pct_rng))
	Threads.@threads for i = 1:length(r_repo_pct_rng)
		conversions[i] = conversion_am(treeType, S, K, r_csa_mm_pct, r_repo_pct_rng[i], vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
	end

	global p_repo = scatter(r_repo_pct_rng, conversions, label = "conv price as function of repo rate", xlabel = "r_repo_pct", ylabel = "conversion px", color=:blue)
	plot!(r_repo_pct_rng, conversions, label = "", color=:blue)
	hline!([mktConversionPx], label = "mkt conv px = $mktConversionPx")
	Print()
	#png(p1,"/Users/alex/Downloads/png1.png")
end

# ╔═╡ 1bd4a082-fcfb-4ae2-b435-d8ea0e926033
if doRepoGraph
	try
		plot(p_repo)
	catch e
		println("the repo graph barfed - no graph available")
	end
else
	Print()
end

# ╔═╡ 7a389e3c-a5d7-4780-b637-122b88d35c4c
if doRiskGraphs
		if sum(divamt) > 0
			title = "Conversion with strike = \$$K, vol = $(vol_pct)%, expiry = $((dt_expiry - dt_trade).value) days, \ndiv = $divamt, exdiv = $(Dates.value.(dt_exdiv .- dt_trade)) days, r_csa = $(r_csa_mm_pct)%, r_repo = $(r_repo_mm_pct)%"
		else
			title = "Conversion with strike = \$$K, vol = $(vol_pct)%, expiry = $((dt_expiry - dt_trade).value) days, \nr_csa = $(r_csa_mm_pct)%, r_repo = $(round(r_repo_mm_pct,digits = 2))%"
		end

		return_range = (100 .+ collect(-30:1:80)) ./ 100.0
		n = length(return_range)
		δ = Vector{Float64}(undef, n)
		Γ = Vector{Float64}(undef, n)
		Threads.@threads for i = 1:n
			s = return_range[i] * S
			putVal, putDelta, putGamma = option_am(treeType, "put", s, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
			callVal, callDelta, callGamma = option_am(treeType, "call", s, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
			δ[i] = putDelta - callDelta + 1
			Γ[i] = putGamma - callGamma
		end
		
		p1a = plot(
			S .* return_range,
			δ,
			lw = 4,
			label = "",
			ylabel = "Delta",
			#xlabel = "Price",
			#xguidefontsize = 8,
			title = title,
			titlefontsize = 10,
		)
		vline!(p1a, [S], legend = false)
		p1b = plot(
			S .* return_range,
			Γ,
			lw = 4,
			label = "",
			ylabel = "Gamma",
			xlabel = "Price",
			xguidefontsize = 8,
			#title = title,
			#titlefontsize = 10
		)
		vline!(p1b, [S], legend = false)
	Print()
	end

# ╔═╡ c4fe8ad8-bbd9-490f-8509-27a3536cbc2c
if doRiskGraphs
	try
		plot(p1a, p1b, layout = (2, 1))
	catch e
		println("the code barfed - no graph available")
	end
else
	Print()
end

# ╔═╡ 6542d0e2-8e5f-4a06-ac97-4aca938d1594
if doExBoundaryGraphs
	if treeType == "CRR"
		(callpx, X_call) = crr_extra_output("call", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt,2)[1:2]
		#
		(putpx, X_put) = crr_extra_output("put", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt,2)[1:2]
	elseif treeType == "LR"
		if (dt_expiry - dt_trade).value % 2 == 0
			println("Skipping LR valuation because stepsPerDay = 1 gives even number of steps -> unreliable valuation")
			X =[]
		else
			(callpx, X_call) = lr_extra_output("call", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt,2)[1:2]
			#
			(putpx, X_put) = lr_extra_output("put", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt,2)[1:2]
		end
	elseif treeType == "Trinomial"
		(callpx, X_call) = trinomial_extra_output("call", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt,2)[1:2]
		#
		(putpx, X_put) = trinomial_extra_output("put", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt,2)[1:2]
	end

	if sum(divamt) > 0
		dates2 = dt_trade:Day(1):maximum(vcat(dt_divpmt, dt_expiry))
		ndays2 = length(dates2)
		divfvs = zeros(Float64,ndays2)
		# we want at each slice divfv = value of all divs w/ exdiv not in the past as of that slice
		dfs = exp.(-r_csa * (0:ndays2-1) / 360)
		i_exdiv = zeros(Int, length(divamt))
		i_divpmt = zeros(Int, length(divamt))
		for i = 1:length(divamt)
			i_divpmt[i] = findall(dates2 .== dt_divpmt[i])[1]
			i_exdiv[i] = findall(dates2 .== dt_exdiv[i])[1]
			global divfvs[1:i_exdiv[i]-1] .=divfvs[1:i_exdiv[i]-1] + (divamt[i] * dfs[i_divpmt[i]]) ./ dfs[1:i_exdiv[i]-1]
		end
	end

	callMult = 1
	putMult = -1

	##
	myString_call = ""
	if size(X_call)[1] > 0 # i.e. when there is early exercise
		ndays = nrow(X_call)
		#r_csa = log(1.0 + 0.01 * r_csa_mm_pct * ndays / 360) / (ndays / 360)
		#r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / (ndays / 360)

		S0 = S
		X_call.divfvs = zeros(ndays)
		if sum(divamt) > 0
			X_call.divfvs = divfvs[1:ndays]
			S0 = S0 - divfvs[1]
		end
		
    	σ = vol_pct/100 * S / S0
		
		X_call.noSig = S0 .* exp.((r_repo - 0.5 * σ^2) * X_call.sliceNum / 360 ) + X_call.divfvs
		X_call.oneSig = S0 .* exp.((r_repo - 0.5 * σ^2) * X_call.sliceNum / 360 .+ σ * sqrt.(X_call.sliceNum / 360) * 1 * callMult) + X_call.divfvs
		X_call.twoSig = S0 .* exp.((r_repo - 0.5 * σ^2) * X_call.sliceNum / 360 .+ σ * sqrt.(X_call.sliceNum / 360) * 2 * callMult) + X_call.divfvs
		X_call.threeSig = S0 .* exp.((r_repo - 0.5 * σ^2) * X_call.sliceNum / 360 .+ σ * sqrt.(X_call.sliceNum / 360) * 3 * callMult) + X_call.divfvs
		ymin_call = 0.9 * minimum(vcat(X_call.noSig, X_call.threeSig, X_call.sExPrice[isfinite.(X_call.sExPrice)]))
		ymax_call = 1.1 * maximum(vcat(X_call.noSig, X_call.threeSig, X_call.sExPrice[isfinite.(X_call.sExPrice)]))
		gr()
		expStr = Dates.format(dt_expiry,"d-u-yy")
		trdStr = Dates.format(dt_trade,"d-u-yy")	

		indx = findall(isfinite.(X_call.sExPrice))
		if length(indx) > 0
			mindx = minimum(indx)
			myString_call = "   (earliest call exercise boundary at  $(round(X_call.sExPrice[mindx],digits=2)) on   $(X_call.date[mindx]), i.e. $((X_call.date[mindx]-dt_trade).value) days from trade date)"
		end
	end

	##put
	myString_put = ""
	if size(X_put)[1]>0 # in case of no early exercise
		ndays = nrow(X_put)
		#r_csa = log(1.0 + 0.01 * r_csa_mm_pct * ndays / 360) / (ndays / 360)
		#r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / (ndays / 360)

		S0 = S
		X_put.divfvs = zeros(ndays)
	    if sum(divamt) > 0
			X_put.divfvs = divfvs[1:ndays]
			S0 = S0 - divfvs[1]
		end
		
    	σ = vol_pct/100 * S / S0
		
		X_put.noSig = S0 .* exp.((r_repo - 0.5 * σ^2) * X_put.sliceNum / 360 ) + X_put.divfvs
		X_put.oneSig = S0 .* exp.((r_repo - 0.5 * σ^2) * X_put.sliceNum / 360 .+ σ * sqrt.(X_put.sliceNum / 360) * 1 * putMult) + X_put.divfvs
		X_put.twoSig = S0 .* exp.((r_repo - 0.5 * σ^2) * X_put.sliceNum / 360 .+ σ * sqrt.(X_put.sliceNum / 360) * 2 * putMult) + X_put.divfvs
		X_put.threeSig = S0 .* exp.((r_repo - 0.5 * σ^2) * X_put.sliceNum / 360 .+ σ * sqrt.(X_put.sliceNum / 360) * 3 * putMult) + X_put.divfvs
		ymin_put = 0.9 * minimum(vcat(X_put.noSig, X_put.threeSig, X_put.sExPrice[isfinite.(X_put.sExPrice)]))
		ymax_put = 1.1 * maximum(vcat(X_put.noSig, X_put.threeSig, X_put.sExPrice[isfinite.(X_put.sExPrice)]))
		gr()
		expStr = Dates.format(dt_expiry,"d-u-yy")
		trdStr = Dates.format(dt_trade,"d-u-yy")	

		indx = findall(isfinite.(X_put.sExPrice))
		if length(indx) > 0
			mindx = minimum(indx)
			myString_put = "   (earliest put exercise boundary at  $(round(X_put.sExPrice[mindx],digits=2)) on   $(X_put.date[mindx]), i.e. $((X_put.date[mindx]-dt_trade).value) days from trade date)"
		end
	end
	Print()
end

# ╔═╡ 07c8d43e-163d-44dd-ba31-fcf045d178b8
if doExBoundaryGraphs
	gr()
	if size(X_put)[1]>0 # in case of no early exercise
		if false
			p_exput = scatter(X_put.date, X_put.exPrice, ticks=:native, label = "", bg=:lightgray, size = (800,600), ylims = (ymin_put,ymax_put), yguide = "price", yguidefontsize = 12, yguidefontcolor = :black, ytickfontcolor=:black)
		else
			p_exput = scatter(X_put.date, NaN .+ X_put.exPrice, ticks=:native, label = "", bg=:lightgray, size = (800,600), ylims = (ymin_put,ymax_put), yguide = "price", yguidefontsize = 12, yguidefontcolor = :black, ytickfontcolor=:black)
		end
		title!(p_exput, "EE region (light blue) for S = $S, K = $K, $expStr PUT as of $trdStr with σ = $vol_pct%",titlefontsize=9)
		xlabel!(p_exput, "Stock has 16%, 2%, 0.1% prob of being outside respective band")
		#if lowercase(opttype)[1] == 'c'
		#	ribbonGapBelow = zeros(nrow(X_put))
		#	ribbonGapAbove = ymax_call .- X_put.sExPrice
		#else
			ribbonGapBelow2 = X_put.sExPrice .- ymin_put
			ribbonGapAbove2 = zeros(nrow(X_put))	
		#end
		plot!(p_exput, X_put.date, X_put.sExPrice, lw=0, color=:blue, label = "",
			ribbon = (ribbonGapBelow2, ribbonGapAbove2),fillalpha=0.1			
			#ribbon = (X.sExPrice .- ymin,zeros(nrow(X)) ),fillalpha=0.2
		)
		plot!(p_exput, X_put.date, X_put.noSig,label="", color = :black, lw = 2)			
		plot!(p_exput, X_put.date, X_put.oneSig,label="",ls=:dash, color = :darkgray)
		plot!(p_exput, X_put.date, X_put.twoSig,label="",ls=:dash, color = :darkgray)
		plot!(p_exput, X_put.date, X_put.threeSig,label="",ls=:dash, color = :darkgray)
		bar!(twinx(),X_put.date, X_put.cumExProb,ylims=(0,1),color=:red,alpha=0.25,label="", yguide = "cumulative prob of early exercise", yguidefontsize = 12, yguidefontcolor = :red, ytickfontcolor=:red)#yguidefont = font(:red),ytickfont = font(:red, 8),)
	end
	Print()
end

# ╔═╡ 18d602fd-5041-4d08-82f3-74b3efd3b7e5
if doExBoundaryGraphs && earlyExPut
	try
		plot(p_exput)
	catch e
		println("the code barfed - no graph available")
	end
else
	Print()
end

# ╔═╡ fa50e635-ded6-40b3-9a70-b1efd3abaa7c
if doExBoundaryGraphs
	gr()
	if size(X_call)[1]>0 # in case of no early exercise
		if false
			p_excall = scatter(X_call.date, X_call.exPrice, ticks=:native, label = "", bg=:lightgray, size = (800,600), ylims = (ymin_call,ymax_call), yguide = "price", yguidefontsize = 12, yguidefontcolor = :black, ytickfontcolor=:black)
		else
			p_excall = scatter(X_call.date, NaN .+ X_call.exPrice, ticks=:native, label = "", bg=:lightgray, size = (800,600), ylims = (ymin_call,ymax_call), yguide = "price", yguidefontsize = 12, yguidefontcolor = :black, ytickfontcolor=:black)
		end
		title!(p_excall, "EE region (light blue) for S = $S, K = $K, $expStr CALL as of $trdStr with σ = $vol_pct%",titlefontsize=9)
		xlabel!(p_excall, "Stock has 16%, 2%, 0.1% prob of being outside respective band")
		#if lowercase(opttype)[1] == 'c'
			ribbonGapBelow = zeros(nrow(X_call))
			ribbonGapAbove = ymax_call .- X_call.sExPrice
		#else
		#	ribbonGapBelow = X_call.sExPrice .- ymin
		#	ribbonGapAbove = zeros(nrow(X))	
		#end
		plot!(p_excall, X_call.date, X_call.sExPrice, lw=0, color=:blue, label = "",
			ribbon = (ribbonGapBelow, ribbonGapAbove),fillalpha=0.1			
			#ribbon = (X.sExPrice .- ymin,zeros(nrow(X)) ),fillalpha=0.2
		)
		plot!(p_excall, X_call.date, X_call.noSig,label="", color = :black, lw = 2)			
		plot!(p_excall, X_call.date, X_call.oneSig,label="",ls=:dash, color = :darkgray)
		plot!(p_excall, X_call.date, X_call.twoSig,label="",ls=:dash, color = :darkgray)
		plot!(p_excall, X_call.date, X_call.threeSig,label="",ls=:dash, color = :darkgray)
		bar!(twinx(),X_call.date, X_call.cumExProb,ylims=(0,1),color=:red,alpha=0.25,label="", yguide = "cumulative prob of early exercise", yguidefontsize = 12, yguidefontcolor = :red, ytickfontcolor=:red)#yguidefont = font(:red),ytickfont = font(:red, 8),)
		Print()
	end
end

# ╔═╡ e45dcbae-bcb0-4e33-9f0d-94519c47fd04
if doExBoundaryGraphs && earlyExCall
	try
		plot(p_excall)
	catch e
		println("the code barfed - no graph available")
	end
else
	Print()
end

# ╔═╡ 1e467d1d-4b5d-45a1-a6f5-043245a4d728
begin
		if doMtmPnlGraph
			#tic = Sys.time()
			convPxAtTradeDate = mktConversionPx
			if isbad(convPxAtTradeDate)
				convPxAtTradeDate = conversion_Basic(
					S,
					K,
					r_csa_mm_pct,
					r_repo_mm_pct,
					vol_pct,
					dt_exp,
					stepsPerDay,
					dt_exdiv,
					dt_divpmt,
					divamt,
					dt_trade,
					frac_of_trading_day_remaining,
					earlyExPut,
					earlyExCall
				)[1]
			end
	
			rtnVector = collect(LinRange(min_rtn_pct, max_rtn_pct, numPrices))
			synthPxVector = (100 .+ rtnVector) ./ 100 * S
			t_exp = (dt_expiry - dt_trade).value
			dateVector = unique(
				dt_trade +
				Day.(
					round.(
						Int32,
						#collect(LinRange(0, (dt_exp - dt_trade).value, numDates)),
						collect(LinRange(1, (dt_expiry - dt_trade).value, numDates)),
					),
				),
			)
	
			v1 = repeat(dateVector, outer = length(synthPxVector)) # dates
			v2 = [(v1[i] - v1[1]).value for i = 1:length(v1)] # days from trade date
			v3 = repeat(rtnVector, inner = length(dateVector)) # returns
			v4 = repeat(synthPxVector, inner = length(dateVector)) # synthetic pxs
			v5 = copy(v4)
			for i=1:length(divamt)
				local indx = findall(v1 .>= dt_exdiv[i])
				if length(indx) > 0
					v5[indx] = v5[indx] .- divamt[i]
				end
			end
	
			numEval = numPrices * length(dateVector)
	
			pnl_array = zeros(length(dateVector), numPrices)
	
			data = DataFrame(
				date = v1,
				days_out = v2,
				pct_move = v3,
				px_synth = v4,
				px_act = v5,
				convPx = zeros(numEval),
				call_ex = (zeros(numEval) .> 1),          #just a way of saying false
				call_ex_days_out = zeros(numEval) .+ NaN, #just a way of saying NaN
				call_ex_px = zeros(numEval) .+ NaN,
				unex_put_px = zeros(numEval) .+ NaN,
				put_ex = (zeros(numEval) .> 1),
				put_ex_days_out = zeros(numEval) .+ NaN,
				put_ex_px = zeros(numEval) .+ NaN,
				unex_call_px = zeros(numEval) .+ NaN,
				accruedInterest = zeros(numEval),
				pnl = zeros(numEval),
			)
	
			#=(call, put, exInfo) = conversion(
				 S,
				 K,
				 r_csa,
				 r_repo,
				 vol,
				 dt_exp,
				 steps_per_day,
				 dt_exdiv,
				 dt_divpmt,
				 divamt,
				 dt_trade,
				 true,
				 frac_of_trading_day_remaining,
			 )=#
	
			allDates = exInfo.date#dt_trade .+ Day.(0:(dt_exp-dt_trade).value)
			cb = NaN .+ zeros(length(allDates))
			pb = NaN .+ zeros(length(allDates))
			indx1 = round.(Int32, myMatch(allDates, exInfo.date))
			for i = 1:length(indx1)
				j = indx1[i]
				cb[j] = exInfo.cb[i]
				pb[j] = exInfo.pb[i]
			end
	
			# fill data dataframe with exercise info for linear px evolutions
			t_exdiv_vec = Dates.value.(dt_exdiv .- dt_trade)
			Threads.@threads for i = 1:numEval
				#println(i)
				if data.days_out[i] > 0
					#global px = collect(LinRange(S, data.px_synth[i], 1 + data.days_out[i]))
					#global N = length(px)
					px = collect(LinRange(S, data.px_synth[i], 1 + data.days_out[i]))
					N = length(px)
					# adjust p for divs
					for j=1:length(t_exdiv_vec)
						if N >= 1 + t_exdiv_vec[j]
							px[1+t_exdiv_vec[j]:N] = px[1+t_exdiv_vec[j]:N] .- divamt[j]
						end
					end
					#global pindx = findall(px .<= pb[1:N])
					pindx = findall(px .<= pb[1:N])
					if length(pindx) > 0
						pindx = minimum(pindx)
						data.put_ex[i] = true
						data.put_ex_px[i] = px[pindx]
						data.put_ex_days_out[i] = pindx - 1
						if useAmericanPricerForUnexOpts
							data.unex_call_px[i] = option_am(
								"call",
								data.put_ex_px[i],
								K,
								r_csa_mm_pct,
								r_repo_mm_pct,
								vol_pct,
								dt_exp,
								5,
								dt_exdiv,
								dt_divpmt,
								divamt,
								dt_trade + Day(data.put_ex_days_out[i]),
								frac_of_trading_day_remaining,
							)
						else
							data.unex_call_px[i] = BS(
								"call",
								data.put_ex_px[i],
								K,
								r_csa_mm_pct,
								r_repo_mm_pct,
								vol_pct,
								dt_exp,
								dt_exdiv,
								dt_divpmt,
								divamt,
								dt_trade + Day(data.put_ex_days_out[i]),
							)
						end
					end
					#global cindx = findall(px .>= cb[1:N])
					cindx = findall(px .>= cb[1:N])
					if length(cindx) > 0
						cindx = minimum(cindx)
						data.call_ex[i] = true
						data.call_ex_px[i] = px[cindx]
						data.call_ex_days_out[i] = cindx - 1
						if useAmericanPricerForUnexOpts
							data.unex_put_px[i] = option_am(
								"put",
								data.call_ex_px[i],
								K,
								r_csa_mm_pct,
								r_repo_mm_pct,
								vol_pct,
								dt_exp,
								5,
								dt_exdiv,
								dt_divpmt,
								divamt,
								dt_trade + Day(data.call_ex_days_out[i]),
								frac_of_trading_day_remaining,
							)
						else
							data.unex_put_px[i] = BS(
								"put",
								data.call_ex_px[i],
								K,
								r_csa_mm_pct,
								r_repo_mm_pct,
								vol_pct,
								dt_exp,
								dt_exdiv,
								dt_divpmt,
								divamt,
								dt_trade + Day(data.call_ex_days_out[i]),
							)
						end
					end
				end
			end
	
			Threads.@threads for i = 1:numEval
				#println(i)
				a = findfirst(dateVector .== data.date[i])
				b = findfirst(rtnVector .== data.pct_move[i])
	
				if data.date[i] == dt_exp
					data.convPx[i] = 0
					S_avg = (S + data.px_synth[i]) / 2
					data.accruedInterest[i] =
						-(S_avg * r_repo_mm_pct/100 + (K - S_avg + convPxAtTradeDate) * r_csa_mm_pct/100) *
						data.days_out[i] / 360
					data.pnl[i] =
						-convPxAtTradeDate +
						+data.accruedInterest[i] +
						sum(divamt .* (data.days_out[i] .>= t_exdiv_vec))
				elseif data.date[i] < dt_exp
					data.convPx[i] = conversion_Basic(
						data.px_act[i],
						K,
						r_csa_mm_pct,
						r_repo_mm_pct,
						vol_pct,
						dt_exp,
						stepsPerDay,
						dt_exdiv,
						dt_divpmt,
						divamt,
						data.date[i],
					)[1]
					## assumes interest calc'd on avg px
					if !(data.put_ex[i] | data.call_ex[i])
						# pay K + convPx for conversion bundle
						# repo out the stock and use proceeds of S to defray borrowing of K + convPx
						# therefore pay r_csa on K-S + convPx and pay r_repo on S
						S_avg = (data.px_synth[i] + S) / 2
						data.accruedInterest[i] =
							-(
								S_avg * r_repo_mm_pct/100 +
								(K - S_avg + convPxAtTradeDate) * r_csa_mm_pct/100
							) * data.days_out[i] / 360
						data.pnl[i] =
							(data.convPx[i] - convPxAtTradeDate) +
							data.accruedInterest[i] +
							sum(divamt .* (data.days_out[i] .>= t_exdiv_vec))
					else
						if data.put_ex[i]
							t_early = data.put_ex_days_out[i]
							S_ex = data.put_ex_px[i]
						else
							t_early = data.call_ex_days_out[i]
							S_ex = data.call_ex_px[i]
						end
						S_avg = (S + S_ex) / 2
						data.accruedInterest[i] =
							-(
								S_avg * r_repo_mm_pct/100 +
								(K - S_avg + convPxAtTradeDate) * r_csa_mm_pct/100
							) * t_early / 360
						# we ignore any residual value in the unexercised option
						data.pnl[i] =
							-convPxAtTradeDate +
							+data.accruedInterest[i] +
							sum(divamt .* (t_early .>= t_exdiv_vec))
						if data.put_ex[i] &
						   includeUnexOptVals &
						   (data.unex_call_px[i] > 0)
							data.pnl[i] = data.pnl[i] - data.unex_call_px[i]
						end
						if data.call_ex[i] &
						   includeUnexOptVals &
						   (data.unex_put_px[i] > 0)
							data.pnl[i] = data.pnl[i] + data.unex_put_px[i]
						end
					end
				end
				pnl_array[a, b] = data.pnl[i]
			end
	
	
			T = zeros(length(dateVector))
			for i = 1:length(T)
				T[i] = (dateVector[i] - dateVector[1]).value
			end
			up1 =
					100 .* (exp.((r_repo - 0.5 * vol^2) .* T / 365 .+ vol .* sqrt.(T ./ 365) * 1) .- 1)
			up2 =
				100 .* (exp.((r_repo - 0.5 * vol^2) .* T / 365 .+ vol .* sqrt.(T ./ 365) * 2) .- 1)
			up3 =
				100 .* (exp.((r_repo - 0.5 * vol^2) .* T / 365 .+ vol .* sqrt.(T ./ 365) * 3) .- 1)
			down1 =
				100 .* (exp.((r_repo - 0.5 * vol^2) .* T / 365 .+ vol .* sqrt.(T ./ 365) * -1) .- 1)
			down2 =
				100 .* (exp.((r_repo - 0.5 * vol^2) .* T / 365 .+ vol .* sqrt.(T ./ 365) * -2) .- 1)
			down3 =
				100 .* (exp.((r_repo - 0.5 * vol^2) .* T / 365 .+ vol .* sqrt.(T ./ 365) * -3) .- 1)
			centerLine = zeros(length(T))
	
	
			if sum(divamt) > 0
				xlabel6 = "expiry = $dt_exp, K = $K, div = $divamt, exdiv_dt = $dt_exdiv, r_csa = $r_csa_mm_pct%, rebate = $r_repo_mm_pct%"
			else
				xlabel6 = "expiry = $dt_exp, K = $K, r_csa = $r_csa_mm_pct%, rebate = $r_repo_mm_pct%"
			end
	
			amax = maximum(abs.(pnl_array))
			p_contour = contour(
				dateVector,
				rtnVector,
				transpose(pnl_array),
				fill = (true, cgrad(:RdYlBu_11)),
				clims = (-amax, amax),
			)
			plot!(dateVector, up1, color = :black, legend = false, linestyle = :dot)
			plot!(dateVector, up2, color = :black, legend = false, linestyle = :dot)
			plot!(dateVector, up3, color = :black, legend = false, linestyle = :dot)
			plot!(dateVector, down1, color = :black, legend = false, linestyle = :dot)
			plot!(dateVector, down2, color = :black, legend = false, linestyle = :dot)
			plot!(dateVector, down3, color = :black, legend = false, linestyle = :dot)
			plot!(
				dateVector,
				centerLine,
				color = :black,
				legend = false,
				linestyle = :solid,
			)
			#xlabel!("Date")
			ylabel!("Percentage Move")
			xlabel!(xlabel6, guidefontsize = 7)
			title!(
				"Conversion MTM PNL Heatmap with 1σ, 2σ, 3σ bands \n(stock has 16%, 2%, 0.1% prob of being outside band)",
				titlefontsize = 10,
			)
			#=
				if do3dGraph
					r1 = round(100 * r_csa, digits = 2)
					r2 = round(100 * r_repo, digits = 2)
					if divamt > 0
						myTitle = "$graphTitle<br>convPx = $(round(convPxAtTradeDate,digits=2))  S = $S  K = $K  expiry = $dt_exp  r_csa = $r1%  rebate = $r2%<br> div = $divamt  exdiv_dt = $dt_exdiv   dt_divpmt = $dt_divpmt"
					else
						myTitle = "$graphTitle<br>convPx = $(round(convPxAtTradeDate,digits=2))  S = $S  K = $K  expiry = $dt_exp  r_csa = $r1%  rebate = $r2%"
					end
					plotly()
					p7 = plot(dateVector, rtnVector, transpose(pnl_array), size = (1000,700),st = :surface,title = myTitle)
				end
				Print()
				=#
				#toc = Sys.time()
				#println("done with the mtmContourGraph section in $(round(toc-tic,digits =1)) seconds")
		end
		Print()
	
end

# ╔═╡ 3049b048-2929-401d-8687-1e80ca2b8f48
if doMtmPnlGraph
	try
		plot(p_contour)
	catch e
		Print("the code barfed - no graph available")
	end
else
	Print()
end

# ╔═╡ b6aaa574-7002-42dd-b195-765bb88896c6
begin
	#=
	  if Sys.isapple()
		myPath = "/Users/alex/Library/CloudStorage/Dropbox/Code/Julia/"
     elseif Sys.islinux()
		 myPath = "/home/alex/Dropbox/Code/Julia/"
     end
	using Revise, CSV, DataFrames, Dates, Distributions, Format, Grep, HypertextLiteral, Ipopt, JuMP, LinearAlgebra, Loess, NaNStatistics, Optim, Plots, PlutoUI, PrettyTables, Random, Revise, XLSX
	include(myPath * "myFuncs.jl")
	Print()
	=#
end

# ╔═╡ Cell order:
# ╟─03142ff5-32f4-44a7-a11f-38d55ce06712
# ╟─66f07370-1190-4f2a-a509-fb0bdc72653d
# ╟─bd181597-a3dc-4a15-91f0-9611a330c8db
# ╟─9f048c1e-8e5c-4e78-a587-3e8566fdc098
# ╟─707e1bdc-e84b-4d5c-b44a-a6a61b877e85
# ╟─877df67b-6950-4a5c-881c-b252c77e2707
# ╟─d48514f9-8e54-4666-a0ae-4a5a62daae2f
# ╟─0c482b6d-cd42-4036-8cea-7582c774eef3
# ╟─3046063a-0396-4e8c-a602-b9f8f17295d4
# ╟─3049b048-2929-401d-8687-1e80ca2b8f48
# ╟─18d602fd-5041-4d08-82f3-74b3efd3b7e5
# ╟─e45dcbae-bcb0-4e33-9f0d-94519c47fd04
# ╟─c4fe8ad8-bbd9-490f-8509-27a3536cbc2c
# ╟─1bd4a082-fcfb-4ae2-b435-d8ea0e926033
# ╟─6413cc25-ac66-4989-8307-312a35243d1f
# ╟─97498b0f-b148-41e7-b16e-abb6447f0b61
# ╟─7a389e3c-a5d7-4780-b637-122b88d35c4c
# ╟─6542d0e2-8e5f-4a06-ac97-4aca938d1594
# ╟─07c8d43e-163d-44dd-ba31-fcf045d178b8
# ╟─fa50e635-ded6-40b3-9a70-b1efd3abaa7c
# ╟─1e467d1d-4b5d-45a1-a6f5-043245a4d728
# ╟─bde7c9cc-a124-419b-ab02-6c2273e87a26
# ╟─f5c2fd63-5fc6-40bd-9a7f-9731ea5636b8
# ╟─39493c2e-0870-417b-94bc-3658c1c18d34
# ╟─c1003ee9-c023-41ea-a4ef-b1493287afa3
# ╟─b6aaa574-7002-42dd-b195-765bb88896c6
# ╟─a784bc70-4b47-4154-a2c2-56b676ed91e8
# ╟─86d497c4-109a-4642-996d-50fb9bb6c6d0
