
#using Revise, CSV, DataFrames, Dates, Distributions, Format, Grep, HypertextLiteral, Ipopt, JuMP, LinearAlgebra, Loess, NaNStatistics, Optim, Plots, PrettyTables, Random, XLSX

# holiday calendars from https://www.sifma.org/resources/general/us-holiday-archive/
################ Date and calendar functions
us_trading_holidays_2015 = [
    Date(2015, 1, 1),
    Date(2015, 1, 19),
    Date(2015, 4, 3),
    Date(2015, 5, 25),
    Date(2015, 7, 3),
    Date(2015, 9, 7),
    Date(2015, 10, 12),
    Date(2015, 11, 11),
    Date(2015, 11, 26),
    Date(2015, 12, 25),
]
us_trading_holidays_2016 = [
    Date(2016, 1, 1),
    Date(2016, 1, 18),
    Date(2016, 3, 25),
    Date(2016, 5, 30),
    Date(2016, 7, 4),
    Date(2016, 9, 5),
    Date(2016, 10, 10),
    Date(2016, 11, 11),
    Date(2016, 11, 24),
    Date(2016, 12, 26),
]
us_trading_holidays_2017 = [
    Date(2017, 1, 2),
    Date(2017, 1, 16),
    Date(2017, 2, 20),
    Date(2017, 4, 14),
    Date(2017, 5, 29),
    Date(2017, 7, 4),
    Date(2017, 9, 4),
    Date(2017, 10, 9),
    Date(2017, 11, 23),
    Date(2017, 12, 25),
]
us_trading_holidays_2018 = [
    Date(2018, 1, 1),
    Date(2018, 1, 15),
    Date(2018, 2, 19),
    Date(2018, 3, 30),
    Date(2018, 5, 28),
    Date(2018, 7, 4),
    Date(2018, 9, 3),
    Date(2018, 10, 8),
    Date(2018, 11, 12),
    Date(2018, 11, 22),
    Date(2018, 12, 25),
]
us_trading_holidays_2019 = [
    Date(2019, 1, 1),
    Date(2019, 1, 21),
    Date(2019, 2, 18),
    Date(2019, 4, 19),
    Date(2019, 5, 27),
    Date(2019, 7, 4),
    Date(2019, 9, 2),
    Date(2019, 10, 14),
    Date(2019, 11, 11),
    Date(2019, 11, 28),
    Date(2019, 12, 25),
]
us_trading_holidays_2020 = [
    Date(2020, 1, 1),
    Date(2020, 1, 20),
    Date(2020, 2, 17),
    Date(2020, 4, 10),
    Date(2020, 5, 25),
    Date(2020, 7, 3),
    Date(2020, 9, 7),
    Date(2020, 10, 12),
    Date(2020, 11, 11),
    Date(2020, 11, 26),
    Date(2020, 12, 25),
]

us_trading_holidays_2021 = [
    Date(2021, 1, 1),
    Date(2021, 1, 18),
    Date(2021, 2, 15),
    Date(2021, 4, 2),
    Date(2021, 5, 31),
    Date(2021, 7, 5),
    Date(2021, 9, 6),
    Date(2021, 10, 11),
    Date(2021, 11, 11),
    Date(2021, 11, 25),
    Date(2021, 12, 24),
]
us_trading_holidays_2022 = [
    Date(2022, 1, 17),
    Date(2022, 2, 21),
    Date(2022, 4, 15),
    Date(2022, 5, 30),
    Date(2022, 6, 20),
    Date(2022, 7, 4),
    Date(2022, 9, 5),
    Date(2022, 11, 24),
    Date(2022, 12, 26),
]
us_trading_holidays_2023 = [
    Date(2023, 1, 2),
    Date(2023, 1, 16),
    Date(2023, 2, 20),
    Date(2023, 4, 7),
    Date(2023, 5, 29),
    Date(2023, 6, 19),
    Date(2023, 7, 4),
    Date(2023, 9, 4),
    Date(2023, 11, 23),
    Date(2023, 12, 25),
]
us_trading_holidays_2024 = [
    Date(2024, 1, 1),
    Date(2024, 1, 15),
    Date(2024, 2, 19),
    Date(2024, 3, 29),
    Date(2024, 5, 27),
    Date(2024, 6, 19),
    Date(2024, 7, 4),
    Date(2024, 9, 2),
    Date(2024, 11, 28),
    Date(2024, 12, 25),
]
us_trading_holidays_2025 = [
    Date(2025, 1, 1),
    Date(2025, 1, 20),
    Date(2025, 2, 17),
    Date(2025, 4, 18),
    Date(2025, 5, 26),
    Date(2025, 6, 19),
    Date(2025, 7, 4),
    Date(2025, 9, 1),
    Date(2025, 11, 27),
    Date(2025, 12, 25),
]
us_trading_holidays_2026 = [
    Date(2026, 1, 1),
    Date(2026, 1, 19),
    Date(2026, 2, 16),
    Date(2026, 4, 3),
    Date(2026, 5, 25),
    Date(2026, 6, 19),
    Date(2026, 7, 3),
    Date(2026, 9, 7),
    Date(2026, 11, 26),
    Date(2026, 12, 25),
]
us_trading_holidays = vcat(
    us_trading_holidays_2015,
    us_trading_holidays_2016,
    us_trading_holidays_2017,
    us_trading_holidays_2018,
    us_trading_holidays_2019,
    us_trading_holidays_2020,
    us_trading_holidays_2021,
    us_trading_holidays_2022,
    us_trading_holidays_2023,
    us_trading_holidays_2024,
    us_trading_holidays_2025,
    us_trading_holidays_2026
)
fed_end_of_meeting_dates = [
    Date(2015,1,28)
    Date(2015,3,18)
    Date(2015,4,29)
    Date(2015,6,17)
    Date(2015,7,29)
    Date(2015,9,17)
    Date(2015,10,28)
    Date(2015,12,16)
    Date(2016, 1, 27)
    Date(2016, 3, 16)
    Date(2016, 4, 27)
    Date(2016, 6, 15)
    Date(2016, 7, 27)
    Date(2016, 9, 21)
    Date(2016, 11, 2)
    Date(2016, 12, 14)
    Date(2017, 2, 1)
    Date(2017, 3, 15)
    Date(2017, 5, 3)
    Date(2017, 6, 14)
    Date(2017, 7, 26)
    Date(2017, 9, 20)
    Date(2017, 11, 1)
    Date(2017, 12, 13)
    Date(2018, 1, 31)
    Date(2018, 3, 21)
    Date(2018, 5, 2)
    Date(2018, 6, 13)
    Date(2018, 8, 1)
    Date(2018, 9, 26)
    Date(2018, 11, 8)
    Date(2018, 12, 19)
    Date(2020, 1, 29)
    Date(2020, 3, 3)
    Date(2020, 3, 15)
    Date(2020, 4, 29)
    Date(2020, 6, 10)
    Date(2020, 7, 29)
    Date(2020, 9, 16)
    Date(2020, 11, 5)
    Date(2020, 12, 16)
    Date(2016, 1, 27)
    Date(2016, 3, 17)
    Date(2016, 4, 28)
    Date(2016, 6, 16)
    Date(2016, 7, 28)
    Date(2016, 9, 22)
    Date(2016, 11, 3)
    Date(2016, 12, 15)
    Date(2022, 1, 26)
    Date(2022, 3, 16)
    Date(2022, 5, 4)
    Date(2022, 6, 15)
    Date(2022, 7, 27)
    Date(2022, 9, 21)
    Date(2022, 11, 2)
    Date(2022, 12, 14)
    Date(2023, 2, 1)
    Date(2023, 3, 22)
    Date(2023, 5, 3)
    Date(2023, 6, 14)
    Date(2023, 7, 26)
    Date(2023, 9, 20)
    Date(2023, 11, 1)
    Date(2023, 12, 13)
    Date(2024, 1, 31)
    Date(2024, 3, 20)
    Date(2024, 5, 1)
    Date(2024, 6, 12)
    Date(2024, 7, 31)
    Date(2024, 9, 18)
    Date(2024, 11, 7)
    Date(2024, 12, 18)
    Date(2025, 1, 29)
    Date(2025, 3, 19)
    Date(2025, 5, 7)
    Date(2025, 6, 18)
    Date(2025, 7, 30)
    Date(2025, 9, 17)
    Date(2025, 10, 29)
    Date(2025, 12, 10)
    Date(2026, 1, 28)
    Date(2026, 3, 18)
    Date(2026, 4, 29)
    Date(2026, 6, 17)
    Date(2026, 7, 29)
    Date(2026, 9, 16)
    Date(2026, 10, 28)
    Date(2026, 12, 9)
]

function bd(d, holidays=us_trading_holidays, nonBdDir=1)
    while true
        not_biz_day = in(d, holidays) | (dayofweek(d) > 5)
        if not_biz_day
            d = d + Day(nonBdDir)
        else
            break
        end
    end
    d
end

function workday(d, n, holidays=us_trading_holidays)
    if n > 0
        nonBdDir = 1
    else
        nonBdDir = -1
    end
    d = bd(d, holidays, nonBdDir)
    for i = 1:abs(n)
        d = bd(d + Day(nonBdDir), holidays, nonBdDir)
    end
    d
end

function bizDayAdjustedThirdFriday(d)  #i.e. option expiry of month containing date d
    #rolls back to Thurs if Fri is a holiday
    d = floor(d, Month) # returns 1st day of month
    #dayofweek fn returns 1,2,...,7 for Monday, Tuesday,...,Sunday
    x = dayofweek(d)
    if x > 5
        x = 12 - x
    elseif x <= 5
        x = 5 - x
    end
    x = x + 14
    return workday(d + Day(x), 0)
end

################ Fitting and utility functions

function ls_fit(x, y, xhat, order)
    # does least squares fit of y on x and returns
    # (xhat, yhat = f(xhat), tuple of coeffs) where f(x) comes from the
    # least squares fit of y on x
    iter = Iterators.flatten((0, 1:order))
    X = zeros(length(x), length(collect(iter)))
    for (indx, value) in Iterators.enumerate(iter)
        X[:, indx] = x .^ value
    end
    #B = y
    ls_coeffs = try
        X \ y
    catch e
        [NaN, NaN]
    end
    yhat = zeros(length(xhat))
    for (indx, value) in Iterators.enumerate(iter)
        yhat = yhat .+ ls_coeffs[indx] .* xhat .^ value
    end
    return (xhat, yhat, ls_coeffs)
end

function ls_geomean(x, y) #as per the Tofallis paper 
    X = hcat(ones(length(x)), x)
    Y = hcat(ones(length(y)), y)
    (b1, m1) = X\y #coeffs of regressing y on x, y = m1 * x + b + noise
    (bint, mint) = Y\x #coeffs of regressing x on y, i.e. x = mint * y + bint + noise  <--> y = (1/mint) * x + (-bind/mint) + noise
    (b2, m2) = (-bint/mint, 1/mint)
    m = sqrt(m1 * m2)
    b = sum(y - m * x) / length(x)
    ξ = y .- (m*x .+ b)
    return(b, m, ξ)
end

function ls_fit2(x, y, xhat, order)
    #uses only even powers of x plus linear to fit

    # does least squares fit of y on x and returns
    # (xhat, yhat = f(xhat)) where f(x) comes from the
    # least squares fit of y on x
    if isodd(order)
        order = order - 1
    end
    pows = sort(vcat(1, collect(0:2:order)))

    A = zeros(length(x), length(pows))
    for i in eachindex(pows)
        A[:, i] = x .^ pows[i]
    end
    B = y
    ls_coeffs = A \ B
    yhat = zeros(length(xhat))
    for i in eachindex(pows)
        yhat = yhat .+ ls_coeffs[i] .* xhat .^ pows[i]
    end
    return (xhat, yhat, ls_coeffs)
end

################ Yield Curve functions

function updateRatesHist()
    if Sys.isapple()
        path = "/Users/alex/Library/CloudStorage/Dropbox/Code/Julia/data/on_rate_hists/"
    elseif Sys.islinux()
        path = "/home/alex/Dropbox/Code/Julia/data/on_rate_hists/"
    end
    # assumes Search.xlsx has been created from https://www.newyorkfed.org/markets/reference-rates/obfr
    newData = DataFrame(XLSX.readtable(path * "Search.xlsx", "Results"))[:, 1:3]
    fnam1 = path * "sofr_obfr_ff_hist.csv"
    if filesize(fnam1) > 0
        existingData = DataFrame(CSV.File(fnam1))[:, 1:3]
        allData = unique(vcat(newData, existingData))
    else
        allData = newData
    end

    CSV.write(fnam1, allData)

    fnam2 = path * "rateHist.csv"

    rateHist = copy(allData)
    #rateHist = DataFrame(CSV.File(fnam1))
    #rateHist = rateHist[:, 1:3]
    rename!(rateHist, [:date, :variable, :value])
    rateHist = unique(rateHist)
    rateHist = unstack(rateHist)
    rateHist.date = Date.(rateHist.date, dateformat"m/d/y")
    rateHist = sort(rateHist, order(:date, rev=false))
    rename!(rateHist, [:date, :ff, :obfr, :sofr])
    # now create total return series which we'll need for futures inside their accrual periods
    rateHist.sofrTotRet = coalesce.(rateHist.sofr,0) .+ NaN
    rateHist.obfrTotRet = coalesce.(rateHist.obfr,0) .+ NaN
    sofrStartIndx = findfirst(coalesce.(rateHist.sofr,-3) .> 0)
    obfrStartIndx = findfirst(coalesce.(rateHist.obfr, -3) .> 0)
    rateHist.sofrTotRet[sofrStartIndx] = 1
    rateHist.obfrTotRet[obfrStartIndx] = 1
    for r = 2:nrow(rateHist)
        if r > sofrStartIndx && !ismissing(rateHist.sofr[r])
            if !ismissing(rateHist.sofr[r-1])
                rateHist.sofrTotRet[r] = rateHist.sofrTotRet[r-1] * (1.0 + rateHist.sofr[r-1] / 100 * (rateHist.date[r] - rateHist.date[r-1]).value / 360)
            else #the correctness of this relies upon the fact that there aren't 2 sequential missing values
                rateHist.sofrTotRet[r] = rateHist.sofrTotRet[r-2] * (1.0 + rateHist.sofr[r-2] / 100 * (rateHist.date[r] - rateHist.date[r-2]).value / 360)
            end
        end
        if r > obfrStartIndx && !ismissing(rateHist.obfr[r])
            if !ismissing(rateHist.obfr[r-1])
                rateHist.obfrTotRet[r] = rateHist.obfrTotRet[r-1] * (1.0 + rateHist.obfr[r-1] / 100 * (rateHist.date[r] - rateHist.date[r-1]).value / 360)
            else #the correctness of this relies upon the fact that there aren't 2 sequential missing values
                rateHist.obfrTotRet[r] = rateHist.obfrTotRet[r-2] * (1.0 + rateHist.obfr[r-2] / 100 * (rateHist.date[r] - rateHist.date[r-2]).value / 360)
            end
        end
    end
    CSV.write(fnam2, rateHist)
end

function updateFutsHist(eliminateZeroVolumeDays=false)
    
    if Sys.isapple()
        path = "/Users/alex/Library/CloudStorage/Dropbox/Code/Julia/data/"
    elseif Sys.islinux()
        path = "/home/alex/Dropbox/Code/Julia/data/"
    end

    monthCodes = DataFrame(code=["F", "G", "H", "J", "K", "M", "N", "Q", "U", "V", "X", "Z"], month=January:December)

    function thirdwednesday(d::Date)
        d = floor(d, Month) # returns 1st day of month
        #dayofweek fn returns 1,2,...,7 for Monday, Tuesday,...,Sunday
        x = dayofweek(d)
        if 3 - x < 0
            x = 10 - x
        else
            x = 3 - x
        end
        x = x + 14
        return d + Day(x)
    end

    # firstDay will only work for futures after the year 2000
    function firstDay(ticker::String)
        # if 3rd Wednesday or 1st of month is non-biz day then returns prev biz day
        month = monthCodes.month[findall(monthCodes.code .== uppercase(string(ticker[3])))[1]]
        #year = Year(today()).value - Year(today()).value % 10 + parse(Int64, ticker[5]) # only worked for futs in the 2020s!
        year = 2000 + 10 * parse(Int64, ticker[4]) + parse(Int64, ticker[5])
        if ticker[1:2] == "sq"
            startDate = thirdwednesday(Date(year, month, 1))
        elseif ticker[1:2] ∈ ["sl", "zq"]
            startDate = Date(year, month, 1)
        end
        return startDate
    end

    # lasttDay will only work for futures after the year 2000
    function lastDay(ticker::String)
        month = monthCodes.month[findall(monthCodes.code .== uppercase(string(ticker[3])))[1]]
        #year = Year(today()).value - Year(today()).value % 10 + parse(Int64, ticker[5]) # only worked for futs in the 2020s!
        year = 2000 + 10 * parse(Int64, ticker[4]) + parse(Int64, ticker[5])
        if ticker[1:2] == "sq"
            if month < 10
                endDate = thirdwednesday(Date(year, month + 3, 1)) - Day(1)
            else
                endDate = thirdwednesday(Date(year + 1, month - 9, 1)) - Day(1)
            end
        elseif ticker[1:2] ∈ ["sl", "zq"]
            d = Date(year, month, 1)
            endDate = lastdayofmonth(d)
        end
        return endDate
    end

    filepath = path * "futsHist.csv"
    global futsHist = DataFrame(date=String[], last=Float64[], volume=Int64[], openInt=Int64[], bct=String[])#, firstRateDate=Date[], lastRateDate=Date[])

    for dirName ∈ path .* ["sofr_3m_futs/", "sofr_1m_futs/", "ff_futs/"]
        fnams = grep("_daily_historical-data", readdir(dirName))
        for fnam ∈ fnams
            bct = fnam[1:5]
            tmp = DataFrame(CSV.File(dirName * fnam, silencewarnings=true))
            tmp = tmp[1:nrow(tmp)-1, :]
            tmp.bct = string.(1:nrow(tmp))
            tmp.bct .= bct
            tmp = tmp[:, [1, 5, 8, 9, 10]]
            rename!(tmp, [:date, :last, :volume, :openInt, :bct])
            futsHist = vcat(futsHist, tmp)
        end
    end

    function convert_to_dates(date_strings::AbstractVector)
        # First ensure we have Vector{String}
        str_vector = Vector{String}(date_strings)

        dates = similar(str_vector, Date)
        for (i, ds) in enumerate(str_vector)
            if occursin("/", ds)
                # Handle mm/dd/yyyy format
                dates[i] = Date(ds, dateformat"mm/dd/yyyy")
            else
                # Handle yyyy-mm-dd format
                dates[i] = Date(ds, dateformat"yyyy-mm-dd")
            end
        end
        return dates
    end

    futsHist.date .= convert_to_dates(futsHist.date)
    sort!(futsHist,:date)

    
    futsHist.last = Float64.(futsHist.last)
    futsHist.volume = Int64.(futsHist.volume)
    futsHist.openInt = Int64.(futsHist.openInt)
    futsHist.firstRateDate = copy(futsHist.date)
    futsHist.lastRateDate = copy(futsHist.date)
    for r = 1:nrow(futsHist)
        futsHist.firstRateDate[r] = firstDay(futsHist.bct[r])
        futsHist.lastRateDate[r] = lastDay(futsHist.bct[r])
    end
    if eliminateZeroVolumeDays
        futsHist = futsHist[futsHist.volume .>0, :]
    end
    # because we're going to interleave 1m and 3m sofr futs, let's give them a common date range
    indx_sq = findall(contains.(futsHist.bct, "sq") .== 1)
    indx_sl = findall(contains.(futsHist.bct, "sl") .== 1)
    min_date_sq = minimum(futsHist.date[indx_sq])
    max_date_sq = maximum(futsHist.date[indx_sq])
    min_date_sl = minimum(futsHist.date[indx_sl])
    max_date_sl = maximum(futsHist.date[indx_sl])
    min_date = max(min_date_sq, min_date_sl)
    max_date = min(max_date_sq, max_date_sl)
    rows2delete1 = findall(contains.(futsHist.bct, "sq") .== 1 .&& futsHist.date .< min_date)
    rows2delete2 = findall(contains.(futsHist.bct, "sq") .== 1 .&& futsHist.date .> max_date)
    rows2delete3 = findall(contains.(futsHist.bct, "sl") .== 1 .&& futsHist.date .< min_date)
    rows2delete4 = findall(contains.(futsHist.bct, "sl") .== 1 .&& futsHist.date .> max_date)
    delete!(futsHist, sort(unique(union(rows2delete1, rows2delete2, rows2delete3, rows2delete4))))
    sort!(futsHist, [:date, :lastRateDate])

    #if false
    ## change this part so it's not dependent on EST
    lastFullTradingDay = bd(today(), us_trading_holidays, -1)
    c1 = Hour(now()).value >= 10
    c2 = Hour(now()).value == 9 && Minute(now()).value > 30
    c3 = Hour(now()).value <= 15
    if (c1 || c2) && c3
        intraDayFlag = true
    else
        intraDayFlag = false
    end
    #futsHist = DataFrame(CSV.File(path * "futsHist.csv"))
    if lastFullTradingDay < today() || intraDayFlag
        futsHist = futsHist[futsHist.date.<today(), :]
    end
    futsHist = futsHist[futsHist.date .<= today(), :] # to remove post close data
    CSV.write(filepath, futsHist)
    #end
end

function curveStripperFromFuts(baseRate, d, dEnd=d + Year(2), futs=missing, d_lastRateSet=missing, r_lastRateSet=missing, r_currRateSet=missing; eliminateZeroVolumeDays=true)

    #=
    		to convexity adjust the futures prices we should raise the price / lower the
    		implied rate by sig^2 T1 T2 / 2 where the futures interval starts at T1 and
    		ends at T2 (sig is vol for the forward rate and has units time^-3/2 because dr ~ sig * dz).  

    		The sign is obvious from the negative convexity implied from margining:  rates go up, fut pxs go down -> borrow money at higher rate to keep up w/ variation margining rates go down, fut pxs go up -> take money out of margin account and invest at lower rates.  this means that by being long the future you're effectively locking in a lower rate than 100 - futpx, hence we reduce the implied rate by the convexity adjustment.

    		the more rates jiggle the more costly this effect also note [rate]= 1/time = 
    		N.B. for SOFR futures that average a series of daily fixings from T1 to T2 we can 
    		approximate things well with sig^2 * Tbar^2 /2 where Tbar = (T1+T2)/2 and sig refers to imp vol for 1day rates at expiry Tbar
    r_repo    		=#
    begin
        if Sys.isapple()
            path = "/Users/alex/Library/CloudStorage/Dropbox/Code/Julia/data/"
        elseif Sys.islinux()
            path = "/home/alex/Dropbox/Code/Julia/data/"
        end
    end

    #=
    		we will choose lastRateSet and d_lastRateSet to reference d rather than workday(d,-1).  this is because the rate set only changes on fed dates (to a close approximation) and on a fed date they've announced the rate that's going to print the next day, so we already know it and it's baked into the futures closing prices on d
    		=#

    if ismissing(d_lastRateSet) || ismissing(r_lastRateSet)
        rateHist = DataFrame(CSV.File(path * "/on_rate_hists/rateHist.csv"))
        missingIndx = findall(ismissing.(rateHist[:, lowercase(baseRate)]))
        if length(missingIndx) > 0
            deleteat!(rateHist, missingIndx)
        end
        indx_earlier = findall(rateHist.date .< d)
        indx_missing = findall(ismissing.(rateHist[:, lowercase(baseRate)]))
        indx_earlier = setdiff(indx_earlier, indx_missing)
        i_lastRateSet = maximum(indx_earlier)
        d_lastRateSet = rateHist.date[i_lastRateSet]
        r_lastRateSet = rateHist[i_lastRateSet, lowercase(baseRate)]
    end

    if ismissing(r_currRateSet)
        indx = findall(rateHist.date .>= d)
        if length(indx) > 0
            r_currRateSet = rateHist[minimum(indx), lowercase(baseRate)]
        else
            r_currRateSet = r_lastRateSet
        end
    end

    if ismissing(futs)
        futsHist = DataFrame(CSV.File(path * "futsHist.csv"))
        syms = futsHist.bct
        if lowercase(baseRate) == "sofr"
            indx1 = findall(contains.(syms, "sq") .== 1)
            indx2 = findall(contains.(syms, "sl") .== 1)
            indx = sort(unique(vcat(indx1, indx2)))
        elseif lowercase(baseRate) == "ff"
            indx = findall(contains.(syms, "zq") .== 1)
        end
        futsHist = futsHist[indx, :]
        futs = futsHist[findall(futsHist.date .== d), :]# .&& futsHist.volume .> 100), :] # .&& futsHist.openInt .> 100), :]
        if eliminateZeroVolumeDays
            futs = futs[futs.volume .> 0, :]
        end
        futs = futs[futs.lastRateDate .> d .&& futs.firstRateDate .<= dEnd, :]
        sort!(futs, [:date, :lastRateDate])
        #= we'll build the curve from the first 5 1m sofr futs and then 3m sofr futures on the HMUZ cycle that end after the 5th 1m future ends

        	         *** CME methodology for term SOFR construction is described here:
        			https://www.cmegroup.com/market-data/files/cme-term-sofr-reference-rates-benchmark-methodology.pdf which is, in turn, based on SSRN-id3352598.pdf
        			=#
        if lowercase(baseRate) == "sofr"
            delete!(futs, findall(contains.(futs.bct, "sq") .== 1 .&& Dates.monthname.(futs.lastRateDate) .∉ Ref(["March", "June", "September", "December"])))
        end
    end
    if lowercase(baseRate) == "sofr"
        indx1m = findall(contains.(futs.bct, "sl") .== 1)
        indx1m = indx1m[1:5]
        indx3m = findall(contains.(futs.bct, "sq") .== 1 .&& futs.lastRateDate .> futs.lastRateDate[indx1m[end]])
        futs = futs[sort(vcat(indx1m, indx3m)), :]
    end

    futs.lockedInTRorSum = futs.last .+ NaN

    for r = 1:nrow(futs)
        if futs.firstRateDate[r] < d
            if futs.bct[r][1:2] == "sq"
                # dates ends on last biz day prior to d
                dates = sort(unique(bd.(futs.firstRateDate[r]:Day(1):d-Day(1), Ref(us_trading_holidays), Ref(-1))))
                nd = length(dates)
                prevRateSettings = DataFrame(date=dates, rate=zeros(nd) .+ NaN)
                for i = 1:nd
                    if i == 1
                        j = maximum(findall(rateHist.date .<= dates[1]))
                    else
                        j = findall(rateHist.date .== dates[i])
                    end
                    if length(j) > 0
                        prevRateSettings.rate[i] = rateHist[j[1], lowercase(baseRate)]
                    else # this else clause shouldn't normally be activated since we're looping over biz days only
                        prevRateSettings.rate[i] = prevRateSettings[i-1, 2]
                    end
                end
                futs.lockedInTRorSum[r] = 1
                for i = 1:nd-1
                    futs.lockedInTRorSum[r] = futs.lockedInTRorSum[r] * (1 + 0.01 * prevRateSettings.rate[i] * (dates[i+1] - dates[i]).value / 360)
                end
                futs.lockedInTRorSum[r] = futs.lockedInTRorSum[r] * (1 + 0.01 * prevRateSettings.rate[end] * (d - dates[end]).value / 360)
            else # sl and zq
                dates = futs.firstRateDate[r]:Day(1):d-Day(1)
                nd = length(dates)
                prevRateSettings = DataFrame(date=dates, rate=zeros(nd) .+ NaN)
                for i = 1:nd
                    if i == 1
                        j = maximum(findall(rateHist.date .<= dates[1]))
                    else
                        j = findall(rateHist.date .== dates[i])
                    end
                    if length(j) > 0
                        prevRateSettings.rate[i] = rateHist[j[1], lowercase(baseRate)]
                    else
                        prevRateSettings.rate[i] = prevRateSettings[i-1, 2]
                    end
                end
                futs.lockedInTRorSum[r] = sum(prevRateSettings.rate)
            end
        end
    end

    fedDates = workday.(fed_end_of_meeting_dates, 1)
    fedDates = sort(fedDates[fedDates.>=d.&&fedDates.<=futs.lastRateDate[end]])

    # let's create the intervals within which the o/n rate will be piecewise const

    intervals = DataFrame(startDate=fedDates[1:end-1], endDate=vcat(fedDates[2:end] .- Day(1)))
    indx = findall(futs.lastRateDate .> fedDates[end])
    if length(indx) > 0
        for i ∈ indx
            intervals = vcat(intervals, DataFrame(startDate=intervals.endDate[end] + Day(1), endDate=futs.lastRateDate[i]))
        end
    end
    nVars = nrow(intervals)

    #=
    	    we define 'x' to be the set of piecewise constant o/n rates that we will solve
    		for, and we initilize x to be the average futures rate - just to have something easy and not entirely crazy
    		=#
    xInit = ones(nVars) * (100 - sum(futs.last) / nrow(futs))

    function makeRates(x, rinfo=intervals, ron=r_currRateSet, valDate=d; bizDaysOnly=false)
        rinfo.rate = x
        if valDate < rinfo.startDate[1]
            rinfo = vcat(DataFrame(startDate=valDate, endDate=rinfo.startDate[1] - Day(1), rate=ron), rinfo)
        end

        allDates = valDate:Day(1):rinfo.endDate[end]
        if bizDaysOnly
            allDates = sort(unique(bd.(allDates)))
        end
        nd = length(allDates)
        obj = DataFrame(date=allDates, rate=zeros(nd) .+ NaN)
        for i = 1:nrow(rinfo)
            indx = findall(allDates .>= rinfo.startDate[i] .&& allDates .<= rinfo.endDate[i])
            obj.rate[indx] .= rinfo.rate[i]
        end
        obj
    end

    function calcFutPx(x, symbol, rinfo, futs, r_curr)
        r = findall(string.(futs.bct) .== symbol)[1]
        d = futs.date[1]
        d1 = futs.firstRateDate[r]
        d2 = futs.lastRateDate[r]
        rinfo.rate = x
        if d < rinfo.startDate[1]
            rinfo = vcat(DataFrame(startDate=d, endDate=rinfo.startDate[1] - Day(1), rate=r_curr), rinfo)
        end
        lastRow = minimum(findall(d2 .<= rinfo.endDate))
        firstRow = maximum(findall(d1 .>= rinfo.startDate), init=1)

        if symbol[1:2] == "sq"
            if futs.firstRateDate[r] < d
                tr = futs.lockedInTRorSum[r]
            else
                tr = 1
            end

            for r = firstRow:lastRow
                tr = tr * (1 + 0.01 * rinfo.rate[r] / 360)^(1 + (min(d2, rinfo.endDate[r]) - max(d1, rinfo.startDate[r])).value)
            end
            tr = tr * (1 + 0.01 * rinfo.rate[lastRow] / 360)^((workday(d2, 1) - d2).value - 1)

            rfut = 100 * (tr - 1) * 360 / (workday(d2, 1) - d1).value
            futpx = 100 - rfut
            return futpx
        elseif symbol[1:2] ∈ ["sl", "zq"]
            if futs.firstRateDate[r] < d
                mysum = futs.lockedInTRorSum[r]
            else
                mysum = 0
            end
            # we rely on the fact that only monthly contracts use the average rate to settle
            # therefore, 1) futs.lastDate is the end-of-month and we can just use the "day" function
            # 2) there are at most 2 rows of rinfo that apply

            for r = firstRow:lastRow
                mysum = mysum + rinfo.rate[r] * (1 + (min(d2, rinfo.endDate[r]) - max(d1, rinfo.startDate[r])).value)
            end
            avgrate = mysum / day(d2)
            futpx = 100 - avgrate
        end
        return futpx
    end

    function fmin(x, futs, rinfo, lambda)
        futs = futs[futs.lastRateDate.<=rinfo.endDate[end], :]
        nf = nrow(futs)
        v1 = futs.last
        v2 = zeros(nf) .+ NaN
        for i = 1:nf
            v2[i] = calcFutPx(x, string(futs.bct[i]), rinfo, futs, r_currRateSet)
        end
        1 / nrow(futs) * sum((v1 - v2) .^ 2) + lambda * sum(diff(x) .^ 2)
        #1 / nrow(futs) * sum((v1 - v2) .^ 2) + lambda * sum(diff(x[1:nrow(rinfo)]) .^ 2)
    end


    lambda = 0.001
    res = Optim.optimize(x -> fmin(x, futs, intervals, lambda), xInit, BFGS())#, 	Optim.Options(x_tol=.005))	
    x = Optim.minimizer(res)
    # ycinfo will have rates for all dates but dfs only for biz days
    #makeRates(x, rinfo = intervals, ron=r_currRateSet, valDate=d; bizDaysOnly=false)
    ycinfo = makeRates(x, bizDaysOnly=false)
    ycinfo2 = makeRates(x, bizDaysOnly=true)
    ycinfo2.df = ones(nrow(ycinfo2))
    ycinfo2.df[2:end] = 1 ./ cumprod(1 .+ 0.01 * ycinfo2.rate[1:end-1] .* Dates.value.(diff(ycinfo2.date)) / 360)
    ycinfo.df = zeros(nrow(ycinfo)) .+ NaN
    for r = 1:nrow(ycinfo2)
        #println(ycinfo2.date[r])
        i = findall(ycinfo.date .== ycinfo2.date[r])[1]
        ycinfo.df[i] = ycinfo2.df[r]
    end
    futs.theoLast = futs.last .+ NaN
    for r = 1:nrow(futs)
        futs.theoLast[r] = calcFutPx(x, futs.bct[r], intervals, futs, r_currRateSet)
    end
    futs.bpsError = round.(100 * (futs.theoLast - futs.last), digits=1)
    futs = futs[!, [:date, :bct, :last, :theoLast, :bpsError, :volume, :openInt, :firstRateDate, :lastRateDate, :lockedInTRorSum]]

    intervals.rate = x
    if d < intervals.startDate[1]
        intervals = vcat(DataFrame(startDate=d, endDate=intervals.startDate[1] - Day(1), rate=r_currRateSet), intervals)
    end

    #return ycinfo, futs[!, 1:ncol(futs)-1], intervals
    return ycinfo, futs, intervals
end

function impliedFutPx(symbol, yc)
    (ycinfo, futs, intervals) = yc
    r = findall(string.(futs.bct) .== symbol)[1]
    d = futs.date[1]
    firstRateDate = futs.firstRateDate[r]
    lastRateDate = futs.lastRateDate[r]
    if lowercase(symbol)[1:2] == "sq"
        accrualEndDate = min(ycinfo.date[end], workday(lastRateDate, 1))
        r2 = findall(ycinfo.date .== accrualEndDate)[1]
        if firstRateDate < d
            r1 = 1
            tr = futs.lockedInTRorSum[r]
        else
            r1 = findall(ycinfo.date .== bd(firstRateDate))[1]
            tr = 1
        end
        tr = tr * ycinfo.df[r1] / ycinfo.df[r2]
        rfut = 100 * (tr - 1) * 360 / (ycinfo.date[r2] - firstRateDate).value
        futpx = 100 - rfut
    elseif lowercase(symbol)[1:2] ∈ ["sl", "zq"]
        if firstRateDate < d
            mysum = futs.lockedInTRorSum[r]
        else
            mysum = 0
        end
        indx = findall(ycinfo.date .>= firstRateDate .&& ycinfo.date .<= lastRateDate)
        mysum = mysum + sum(ycinfo.rate[indx])
        avgrate = mysum / day(futs.lastRateDate[r])
        futpx = 100 - avgrate
    end
    return futpx
end

function df(dt, info)
    dt = bd(dt)
    (ycinfo, futs, intervals) = info
    return ycinfo.df[findall(ycinfo.date .== dt)[1]]
end

function fwd_rate(startDate, endDate, info)
    startDate = bd(startDate)
    endDate = bd(endDate)
    df_start = df(startDate, info)
    df_end = df(endDate, info)
    fwd = 100 * (df_start / df_end - 1) * 360 / (endDate - startDate).value
    return fwd
end

function Rmm2Rcc(rmm, t)
    # t in days, rates are in pct
    if t == 0
        return (0)
    elseif rmm <= -36000 / t
        return ("rmm is too negative - drives value negative over t days")
    else
        rcc = 100 * 365 / t * log(1 + 0.01 * rmm * t / 360)
        return (rcc)
    end
end

function Rcc2Rmm(rcc, t)
    # t in days, rates are in pct
    if t == 0
        return (0)
    else
        rmm = 100 * 360 / t * (exp(0.01 * rcc * t / 365) - 1)
        return (rmm)
    end
end

################ European option pricing functions

function BS(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
    # note that divamt, dt_exdiv and dt_divpmt need to be array even if there's only one element
    premium_settlement = 1
    if dt_trade < Date(2024,5,28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end
    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)
    dt_exp_settle = workday(dt_expiry, exercise_settlement)

    days = (dt_expiry - dt_trade).value
    t = days / 365
    r_csa = log(1 + 0.01 * r_csa_mm_pct * days / 360) / t
    r_repo = log(1 + 0.01 * r_repo_mm_pct * days / 360) / t
    σ = vol_pct / 100
    Z_csa = exp(-r_csa * (dt_exp_settle - dt_option_settle).value / 365)
    Z_repo = exp(-r_repo * (dt_exp_settle - dt_stock_settle).value / 365)

    divPV = 0
    for i = 1:length(divamt)
        divPV = divPV + divamt[i] * exp(-r_csa * (dt_divpmt[i] - dt_trade).value / 365)
    end

    S0 = S - divPV
    σ = σ * S / S0
    d1 = (log(S0 / K) + (r_repo + σ^2 / 2) * t) / (σ * √t)
    #d1 = (log(S0 / K) + (r_csa + σ^2 / 2) * t) / (σ * √t)
    d2 = d1 - σ * √t
    dist = Normal()
    if lowercase(optType)[1] == 'c'
        V = Z_csa * (S0 / Z_repo * cdf(dist, d1) - K * cdf(dist, d2))
        #=
        Δ =
        Γ =
        Θ =
        vega = 
        =#
    elseif lowercase(optType)[1] == 'p'
        V = Z_csa * (K * cdf(dist, -d2) - S0 / Z_repo * cdf(dist, -d1))
        #=
        Δ =
        Γ =
        Θ =
        vega = 
        =#
    end
    return V
    #return (V, Δ, Γ, Θ, vega)
end

function implied_vol_BS(optPrice, optType, S, K, r_csa_mm_pct, r_repo_mm_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
    obj(x) = (BS(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, x, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt) - optPrice)^2
    res = Optim.optimize(obj, 0, 1000)
    vol = Optim.minimizer(res)
    vol
end

function conversion_BS(S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
    callPrice = BS("c", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
    putPrice = BS("p", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
    conversion = putPrice - callPrice + S - K
    return (conversion)
end

function conversion_euro(treeType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    if lowercase(treeType) == "crr"
        fn = crr_euro
    elseif lowercase(treeType) == "lr"
        fn = lr_euro
    elseif lowercase(treeType) == "trinomial"
        fn = trinomial_euro
    end

    callPrice = fn("c", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1]
    putPrice = fn("p", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1]
    conversion = putPrice - callPrice + S - K
    return (conversion)
end

function implied_repo_from_conversion_BS(convPrice, S, K, r_csa_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end

    dt_stock_settle = workday(dt_trade, stock_settlement)
    dt_exp_settle = workday(dt_expiry, exercise_settlement)
    t_stock_settle_to_exp_settle = (dt_exp_settle - dt_stock_settle).value

    #= the BS function creates the requisite continuos compounding repo rate as follows:
    days = (dt_expiry - dt_trade).value
    t = days / 365
    r_repo = log(1 + 0.01 * r_repo_mm_pct * days / 360) / t
    we want the argument of log to remain positive so we require
    1 + .01 * r_repo_mm_pct * days / 360 > 0
    r_repo_mm_pct > -360 * 100 / days
    =#

    max_allowable_repo_pct = r_csa_mm_pct + 100
    min_allowable_repo_pct = 1 - 36000 / t_stock_settle_to_exp_settle

    obj(x) = (conversion_BS(S, K, r_csa_mm_pct, x, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt) - convPrice)^2
    res = Optim.optimize(obj, min_allowable_repo_pct, max_allowable_repo_pct)
    if Optim.converged(res)
        return Optim.minimizer(res)
    else
        return NaN
    end
end

function crr_euro(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)

    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end
    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)

    ndays = (dt_expiry - dt_trade).value
    T::Float64 = ndays / 365
    N = ceil(Int, stepsPerDay * ndays) + 2
    Δt = T / (N - 2) # we need to extend the tree by 2 steps to get the right delta, but we don't want to increase the stdev on the terminal slice
    r_csa = log(1.0 + 0.01 * r_csa_mm_pct * ndays / 360) / T
    r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / T

    σ = vol_pct / 100
    divfvs = zeros(Float64, N + 1)
    if sum(divamt) > 0
        # we want at each slice divfv = value of all divs w/ exdiv not in the past as of that slice
        dfs = exp.(-r_csa * floor.(Δt .* collect(-2:N-2)) ./ 365)
        i_exdiv = zeros(Int, length(divamt))
        i_divpmt = zeros(Int, length(divamt))
        for i in eachindex(divamt)
            #consider 2 slices per day: (d-1,d,d,d+1,d+1,d+2,d+2,...)
            # first slice on which d+3 appears is 11th slice:  8 = floor(4 + 2*(3-1))
            i_divpmt[i] = floor(4 + stepsPerDay * ((dt_divpmt[i] - dt_trade).value - 1))
            i_exdiv[i] = floor(4 + stepsPerDay * ((dt_exdiv[i] - dt_trade).value - 1))
            divfvs[1:i_exdiv[i]-1] .+= (divamt[i] * dfs[i_divpmt[i]]) ./ dfs[1:i_exdiv[i]-1]
        end
    end
    # we will discount using R_csa from i_option_settle to end
    # we will grow S0 an R_repo from i_stock_settle to end
    # this means we'll go to the LAST slice on dt_option_settle and dt_stock_settle
    #i_option_settle = 1 + stepsPerDay * (dt_option_settle - dt_trade).value
    #i_stock_settle = 1 + stepsPerDay * (dt_stock_settle - dt_trade).value
    S0 = S - divfvs[1]
    σ = σ * S / S0

    R_csa = exp(r_csa * Δt)
    R_repo = exp(r_repo * Δt)

    U = exp(σ * √Δt)
    D = 1.0 / U
    p = (R_repo - D) / (U - D)
    q = 1.0 - p
    # use pnow,qnow for dt < dt_stock_settle
    pnow = (1.0 - D) / (U - D)
    qnow = 1.0 - pnow
    # use Rnow (i.e. no discounting) for slices with dt < dt_opt_settle
    Rnow = 1.0

    local Z_csa::Float64
    local Z_repo::Float64
    local V::Vector{Float64}
    local p2use::Float64
    local q2use::Float64
    local R2use::Float64
    local px::Float64
    local x::Float64
    local y::Float64
    local tau_csa::Float64
    local tau_repo::Float64
    local dt::Date
    local mult::Float64

    if lowercase(optType)[1] == 'c'
        mult = 1.0
    else
        mult = -1.0
    end

    thetaval1 = thetaval2 = Θ = value = Δ = Γ = 0

    S0 = S0 / (U * D)
    dt = dt_expiry
    dt_exercise_settle = workday(dt, exercise_settlement)
    tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
    tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
    Z_csa = exp(-r_csa * tau_csa)
    Z_repo = exp(-r_repo * tau_repo)
    V = Z_csa * [max(0.0, mult * (S0 * U^i * D^(N - i) / Z_repo - K)) for i = 0:N]
    for n = N-1:-1:0
        dt = dt_trade + Day(ceil((n - 2) / (N - 2) * ndays))
        dt_next = dt_trade + Day(ceil((n + 1 - 2) / (N - 2) * ndays))
        #=
        dt_exercise_settle = workday(dt, exercise_settlement)
        tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
        tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
        Z_csa = exp(-r_csa * tau_csa)
        Z_repo = exp(-r_repo * tau_repo)
        =#
        if dt_next <= dt_stock_settle
            p2use = pnow
        else
            p2use = p
        end
        q2use = 1.0 - p2use
        if dt_next <= dt_option_settle
            R2use = Rnow
        else
            R2use = R_csa
        end
        if dt == bd(dt)
            for i = 0:n
                y = (q2use * V[i+1] + p2use * V[i+2]) / R2use
                #px = S0 * U^i * D^(n - i) + divfvs[n+1]
                #x = Z_csa * mult * (px / Z_repo - K)
                V[i+1] = y#max(x, y)
            end
        else
            for i = 0:n
                V[i+1] = (q2use * V[i+1] + p2use * V[i+2]) / R2use
            end
        end
        if n == 4
            thetaval1 = V[3]
        end
        if n == 2
            Δ = (V[3] - V[1]) / (S0 * U^2 - S0 * D^2)
            Γ = 2 * ((V[3] - V[2]) / (S0 * U^2 - S0 * U * D) - (V[2] - V[1]) / (S0 * U * D - S0 * D^2)) / (S0 * U^2 - S0 * D^2)
            value = V[2]
        end
        if n == 0
            thetaval2 = V[1]
            Θ = (thetaval1 - thetaval2) / (4 * Δt)
            return (value, Δ, Γ, Θ)
        end
    end
end

function lr_euro(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)

    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end
    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)

    ndays = (dt_expiry - dt_trade).value
    T = ndays / 365
    N = ceil(Int, stepsPerDay * ndays) + 2
    if iseven(N)
        N::Int = N + 1
    end
    Δt = T / (N - 2)
    r_csa = log(1.0 + 0.01 * r_csa_mm_pct * ndays / 360) / T
    r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / T

    σ = vol_pct / 100.0
    divfvs = zeros(N + 1)
    if sum(divamt) > 0
        # we want at each slice divfv = value of all divs w/ exdiv not in the past as of that slice
        #dfs = exp.(-r_repo .* floor.(Δt .* collect(0:N)) ./ 365)
        dfs = exp.(-r_csa .* floor.(Δt .* collect(-2:N-2)) ./ 365)

        i_exdiv = zeros(Int, length(divamt))
        i_divpmt = zeros(Int, length(divamt))
        for i = 1:length(divamt)
            i_divpmt[i] = floor(Int, 4 + stepsPerDay * ((dt_divpmt[i] - dt_trade).value - 1))
            i_exdiv[i] = floor(Int, 4 + stepsPerDay * ((dt_exdiv[i] - dt_trade).value - 1))
            divfvs[1:i_exdiv[i]-1] .+= (divamt[i] * dfs[i_divpmt[i]]) ./ dfs[1:i_exdiv[i]-1]
        end
    end
    S0 = S - divfvs[1]
    σ = σ * S / S0

    R_csa = exp(r_csa * Δt)
    R_repo = exp(r_repo * Δt)

    #=
    	U = exp(σ * √Δt)
    	D = 1 / U
    	p = (R_repo - D) / (U - D)
    	q = 1-p
    	=#

    #==#

    S_adj = S0 * exp((r_repo - r_csa) * T)
    d1 = (log(S_adj / K) + T * (r_csa + σ^2 / 2)) / (σ * √T)
    d2 = d1 - σ * √T
    function ppif(z, n)
        1 / 2 + sign(z) / 2 * sqrt(1 - exp(-(n + 1 / 6) * (z / (n + 1 / 3 + 0.1 / (n + 1)))^2))
    end
    p = ppif(d2, N - 2)
    q = 1 - p
    p2 = ppif(d1, N - 2)
    U = exp(r_repo * Δt) * p2 / p
    D = exp(r_repo * Δt) * (1 - p2) / q
    #==#
    # use pnow,qnow for dt < dt_stock_settle
    pnow = (1.0 - D) / (U - D)
    qnow = 1.0 - pnow
    # use Rnow (i.e. no discounting) for slices with dt < dt_opt_settle
    Rnow = 1.0

    local Z_csa::Float64
    local Z_repo::Float64
    local V::Vector{Float64}
    local p2use::Float64
    local q2use::Float64
    local R2use::Float64
    local px::Float64
    local x::Float64
    local y::Float64
    local tau_csa::Float64
    local tau_repo::Float64
    local dt::Date
    local mult::Float64

    if lowercase(optType)[1] == 'c'
        mult = 1.0
    else
        mult = -1.0
    end

    thetaval1 = thetaval2 = Θ = value = Δ = Γ = 0

    S0 = S0 / (U * D)
    dt = dt_expiry
    dt_exercise_settle = workday(dt, exercise_settlement)
    tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
    tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
    Z_csa = exp(-r_csa * tau_csa)
    Z_repo = exp(-r_repo * tau_repo)
    V = Z_csa * [max(0.0, mult * (S0 * U^i * D^(N - i) / Z_repo - K)) for i = 0:N]
    for n = N-1:-1:0
        dt = dt_trade + Day(ceil((n - 2) / (N - 2) * ndays))
        dt_next = dt_trade + Day(ceil((n + 1 - 2) / (N - 2) * ndays))
        #=
        dt_exercise_settle = workday(dt, exercise_settlement)
        tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
        tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
        Z_csa = exp(-r_csa * tau_csa)
        Z_repo = exp(-r_repo * tau_repo)
        =#
        if dt_next <= dt_stock_settle
            p2use = pnow
        else
            p2use = p
        end
        q2use = 1.0 - p2use
        if dt_next <= dt_option_settle
            R2use = Rnow
        else
            R2use = R_csa
        end
        if dt == bd(dt)
            for i = 0:n
                y = (q2use * V[i+1] + p2use * V[i+2]) / R2use
                #px = S0 * U^i * D^(n - i) + divfvs[n+1]
                #x = Z_csa * mult * (px / Z_repo - K)
                V[i+1] = y#max(x, y)
            end
        else
            for i = 0:n
                V[i+1] = (q2use * V[i+1] + p2use * V[i+2]) / R2use
            end
        end
        if n == 4
            thetaval1 = V[3]
        end
        if n == 2
            Δ = (V[3] - V[1]) / (S0 * U^2 - S0 * D^2)
            Γ = 2 * ((V[3] - V[2]) / (S0 * U^2 - S0 * U * D) - (V[2] - V[1]) / (S0 * U * D - S0 * D^2)) / (S0 * U^2 - S0 * D^2)
            value = V[2]
        end
        if n == 0
            thetaval2 = V[1]
            Θ = (thetaval1 - thetaval2) / (4 * Δt)
            return (value, Δ, Γ, Θ)
        end
    end
end

function trinomial_euro(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay, λ=sqrt(2pi) / 2)

    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end
    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)

    ndays = (dt_expiry - dt_trade).value
    T = ndays / 365
    N = ceil(Int, stepsPerDay * ndays) + 1
    Δt = T / (N - 1)
    r_csa = log(1.0 + 0.01 * r_csa_mm_pct * ndays / 360) / T
    r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / T


    divfvs = zeros(Float64, N + 1)
    if sum(divamt) > 0
        # we want at each slice divfv = value of all divs w/ exdiv not in the past as of that slice
        #dfs = exp.(-r_repo * floor.(Δt .* collect(0:N)) ./ 365)
        dfs = exp.(-r_csa * floor.(Δt .* collect(-1:N-1)) ./ 365)
        #dfs = exp.(-r_csa .* Dates.value.(sliceDates .- dt_trade) ./365)
        i_exdiv = zeros(Int, length(divamt))
        i_divpmt = zeros(Int, length(divamt))
        for i in eachindex(divamt)
            i_divpmt[i] = floor(3 + stepsPerDay * ((dt_divpmt[i] - dt_trade).value - 1))
            i_exdiv[i] = floor(3 + stepsPerDay * ((dt_exdiv[i] - dt_trade).value - 1))
            divfvs[1:i_exdiv[i]-1] .+= (divamt[i] * dfs[i_divpmt[i]]) ./ dfs[1:i_exdiv[i]-1]
        end
    end
    # we will discount using R_csa from i_option_settle to end
    # we will grow S0 an R_repo from i_stock_settle to end
    # this means we'll go to the LAST slice on dt_option_settle and dt_stock_settle
    #i_option_settle = 1 + stepsPerDay * (dt_option_settle - dt_trade).value
    #i_stock_settle = 1 + stepsPerDay * (dt_stock_settle - dt_trade).value
    S0 = S - divfvs[1]
    σ = vol_pct / 100
    σ = σ * S / S0

    R_csa = exp(r_csa * Δt)
    R_repo = exp(r_repo * Δt)

    ##
    u = exp(λ * σ * sqrt(Δt))
    m = 1
    d = 1 / u
    M = exp(r_repo * Δt)
    Σ = M^2 * (-1 + exp(σ^2 * Δt))
    pᵤ = (u * (Σ + M^2 - M) - (M - 1)) / ((u - 1) * (u^2 - 1))
    pdown = (u^2 * (Σ + M^2 - M) - u^3 * (M - 1)) / ((u - 1) * (u^2 - 1))
    pₘ = 1 - pᵤ - pdown
    #pₘ = 2 / 3
    #pᵤ = 1 / 6 + (r_repo - σ^2 / 2) * sqrt(Δt / (12 * σ^2))
    # pdown = 1-pm - pu
    #p_u_now = 1 / 6 + (0 - σ^2 / 2) * sqrt(Δt / (12 * σ^2)) #for dt < dt_stock_settle
    p_u_now = (u * (Σ + 1^2 - 1) - (1 - 1)) / ((u - 1) * (u^2 - 1))
    p_d_now = (u^2 * (Σ + 1^2 - 1) - u^3 * (1 - 1)) / ((u - 1) * (u^2 - 1))
    p_m_now = 1 - p_u_now - p_d_now
    ##
    # use Rnow (i.e. no discounting) for slices with dt < dt_opt_settle
    Rnow = 1.0
    local Z_csa::Float64
    local Z_repo::Float64
    local V::Vector{Float64}
    local p2use::Float64
    local q2use::Float64
    local R2use::Float64
    local px::Float64
    local x::Float64
    local y::Float64
    local tau_csa::Float64
    local tau_repo::Float64
    local dt::Date
    local mult::Float64

    if lowercase(optType)[1] == 'c'
        mult = 1.0
    else
        mult = -1.0
    end

    thetaval1 = thetaval2 = Θ = value = Δ = Γ = 0.0

    dt = dt_expiry
    dt_exercise_settle = workday(dt, exercise_settlement)
    tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
    tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
    Z_csa = exp(-r_csa * tau_csa)
    Z_repo = exp(-r_repo * tau_repo)
    V = Z_csa * [max(0.0, mult * (S0 * u^i / Z_repo - K)) for i = -N:N]
    for n = N-1:-1:0
        dt = dt_trade + Day(ceil((n - 1) / (N - 1) * ndays))
        dt_next = dt_trade + Day(ceil((n + 1 - 1) / (N - 1) * ndays))
        #=
        dt_exercise_settle = workday(dt, exercise_settlement)
        tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
        tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
        Z_csa = exp(-r_csa * tau_csa)
        Z_repo = exp(-r_repo * tau_repo)
        =#
        if dt_next <= dt_stock_settle
            p_u = p_u_now
            p_d = p_d_now
            p_m = p_m_now
        else
            p_u = pᵤ
            p_d = pdown
            p_m = pₘ
        end
        if dt_next <= dt_option_settle
            R = Rnow
        else
            R = R_csa
        end
        if dt == bd(dt)
            for i = -n:n
                y = (p_d * V[i+n+1] + p_m * V[i+n+2] + p_u * V[i+n+3]) / R
                #px = S0 * u^i + divfvs[n+1]
                #x = Z_csa * mult * (px / Z_repo - K)
                V[i+n+1] = y#max(x, y)
            end
        else
            for i = -n:n
                V[i+n+1] = (p_d * V[i+n+1] + p_m * V[i+n+2] + p_u * V[i+n+3]) / R
            end
        end
        if n == 2
            thetaval1 = V[3]
        end
        if n == 1
            Δ = (V[3] - V[1]) / (S0 * u - S0 * d)
            Γ = 2 * ((V[3] - V[2]) / (S0 * u - S0 * m) - (V[2] - V[1]) / (S0 * m - S0 * d)) / (S0 * u - S0 * d)
            value = V[2]
        end
        if n == 0
            thetaval2 = V[1]
            Θ = (thetaval1 - thetaval2) / (2 * Δt)
            return (value, Δ, Γ, Θ)
        end
    end
end

function implied_repo_from_conversion_euro(convPrice, treeType, S, K, r_csa_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    #=
     the put/call valuations use r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / T
     so we need to require r_repo_mm_pct > 1 + -36000/ndays
     =#

    T = (dt_expiry - dt_trade).value
    min_repo_pct = 1 - 36000 / T

    global max_allowable_repo_pct = r_csa_mm_pct
    global min_allowable_repo_pct = r_csa_mm_pct

    offset_vec = vcat(30, 0, -1:-1:-3, -10:-10:-100, -200:-100:-900, -1000:-1000:-3000)
    offset_vec = offset_vec[offset_vec.>=(min_repo_pct-r_csa_mm_pct)]
    pvec = zeros(length(offset_vec))
    rvec = r_csa_mm_pct .+ offset_vec
    for i = 1:length(rvec)
        pvec[i] = conversion_euro(treeType, S, K, r_csa_mm_pct, rvec[i], vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
        if (i > 1) && ((pvec[i] - convPrice) * (pvec[i-1] - convPrice) < 0) # i.e. the soln lies in between rvec[i] and rvec[i-1]
            max_allowable_repo_pct = rvec[i-1]
            min_allowable_repo_pct = rvec[i]
            break
        elseif (pvec[i] > convPrice)
            return NaN # pvec is a decreasing fn of r so as we walk r in the negative direction, pvec > mktPrice means we can't cross from below and should exit
        end
    end
    if abs(pvec[end]) > 0
        return NaN # pvec[end] will be non-zero only if no soln was found
    end
    obj(x) = (conversion_euro(treeType, S, K, r_csa_mm_pct, x, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay) - convPrice)^2
    res = Optim.optimize(obj, min_allowable_repo_pct, max_allowable_repo_pct)
    if Optim.converged(res)
        return Optim.minimizer(res)
    end
end

function cev(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_bs_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, beta)
    # note that divamt, dt_exdiv and dt_divpmt need to be array even if there's only one element
    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end
    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)
    dt_exp_settle = workday(dt_expiry, exercise_settlement)

    days = (dt_expiry - dt_trade).value
    T = days / 365
    r_csa = log(1 + 0.01 * r_csa_mm_pct * days / 360) / T
    r_repo = log(1 + 0.01 * r_repo_mm_pct * days / 360) / T
    σ_bs = vol_bs_pct / 100
    Z_csa = exp(-r_csa * (dt_exp_settle - dt_option_settle).value / 365)
    Z_repo = exp(-r_repo * (dt_exp_settle - dt_stock_settle).value / 365)

    divPV = 0
    for i = 1:length(divamt)
        divPV = divPV + divamt[i] * exp(-r_csa * (dt_divpmt[i] - dt_trade).value / 365)
    end

    S0 = S - divPV
    σ_bs = σ_bs * S / S0
    σ = σ_bs * S0^(1 - beta)
    v = σ^2 / (2 * r_repo * (beta - 1)) * (-1 + exp(2 * r_repo * (beta - 1) * T))
    a = ((K * exp(-r_repo))^(2 * (1 - beta))) / (v * (1 - beta)^2)
    b = 1 / (1 - beta)
    c = S0^(2 * (1 - beta)) / (v * (1 - beta)^2)
    q = r_csa - r_repo
    if beta < 1
        if lowercase(optType)[1] == 'c'
            return Z_csa * (S0 / Z_repo * (1 - cdf(NoncentralChisq(b + 2, c), a)) - K * cdf(NoncentralChisq(b, a), c))
        elseif lowercase(optType)[1] == 'p'
            return Z_csa * (K * (1 - cdf(NoncentralChisq(b, a), c)) - S0 / Z_repo * cdf(NoncentralChisq(b + 2, c), a))
        end
    elseif beta > 1
        if lowercase(optType)[1] == 'c'
            return Z_csa * (S0 / Z_repo * (1 - cdf(NoncentralChisq(-b, a), c)) - K * cdf(NoncentralChisq(2 - b, c), a))
        elseif lowercase(optType)[1] == 'p'
            return Z_csa * (K * (1 - cdf(NoncentralChisq(2 - b, c), a)) - S0 / Z_repo * cdf(NoncentralChisq(-b, a), c))
        end
    end
end

################ American option pricing functions

function crr(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)

    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end
    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)

    ndays = (dt_expiry - dt_trade).value
    T::Float64 = ndays / 365
    N = ceil(Int, stepsPerDay * ndays) + 2
    #Δt = T / N
    Δt = T / (N - 2) # we need to extend the tree by 2 steps to get the right delta, but we don't want to increase the stdev on the terminal slice
    r_csa = log(1.0 + 0.01 * r_csa_mm_pct * ndays / 360) / T
    r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / T

    σ = vol_pct / 100.0

    dates = dt_trade:Day(1):dt_expiry

    S0 = S - sum(divamt .* (dt_exdiv .> dt_trade) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt_trade) ./ 360))
    σ = σ * S / S0

    R_csa = exp(r_csa * Δt)
    R_repo = exp(r_repo * Δt)

    U = exp(σ * √Δt)
    D = 1.0 / U
    p = (R_repo - D) / (U - D)
    q = 1.0 - p
    # use pnow,qnow for dt < dt_stock_settle
    pnow = (1.0 - D) / (U - D)
    qnow = 1.0 - pnow
    # use Rnow (i.e. no discounting) for slices with dt < dt_opt_settle
    Rnow = 1.0

    local Z_csa::Float64
    local Z_repo::Float64
    local V::Vector{Float64}
    local p2use::Float64
    local q2use::Float64
    local R2use::Float64
    local px::Float64
    local x::Float64
    local y::Float64
    local tau_csa::Float64
    local tau_repo::Float64
    local dt::Date
    local mult::Float64
    local S0::Float64
    local U::Float64
    local D::Float64
    local divfvs::Vector{Float64}


    if lowercase(optType)[1] == 'c'
        mult = 1.0
    else
        mult = -1.0
    end

    thetaval1 = thetaval2 = Θ = value = Δ = Γ = 0.0
    S0 = S0 / (U * D)
    dt = dt_expiry
    dt_exercise_settle = workday(dt, exercise_settlement)
    tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
    tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
    Z_csa = exp(-r_csa * tau_csa)
    Z_repo = exp(-r_repo * tau_repo)
    V = Z_csa * [max(0.0, mult * (S0 * U^i * D^(N - i) / Z_repo - K)) for i = 0:N]
    for n = N-1:-1:0
        dt = dt_trade + Day(ceil((n - 2) / (N - 2) * ndays))
        dt_next = dt_trade + Day(ceil((n + 1 - 2) / (N - 2) * ndays))
        dt_exercise_settle = workday(dt, exercise_settlement)
        tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
        tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
        Z_csa = exp(-r_csa * tau_csa)
        Z_repo = exp(-r_repo * tau_repo)
        if dt_next <= dt_stock_settle
            p2use = pnow
        else
            p2use = p
        end
        q2use = 1.0 - p2use
        if dt_next <= dt_option_settle
            R2use = Rnow
        else
            R2use = R_csa
        end
        if dt == bd(dt)
            for i = 0:n
                y = (q2use * V[i+1] + p2use * V[i+2]) / R2use
                px = S0 * U^i * D^(n - i) + sum(divamt .* (dt_exdiv .> dt) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt) ./ 360))
                x = Z_csa * mult * (px / Z_repo - K)
                V[i+1] = max(x, y)
            end
        else
            for i = 0:n
                V[i+1] = (q2use * V[i+1] + p2use * V[i+2]) / R2use
            end
        end
        if n == 4
            thetaval1 = V[3]
        end
        if n == 2
            Δ = (V[3] - V[1]) / (S0 * U^2 - S0 * D^2)
            Γ = 2 * ((V[3] - V[2]) / (S0 * U^2 - S0 * U * D) - (V[2] - V[1]) / (S0 * U * D - S0 * D^2)) / (S0 * U^2 - S0 * D^2)
            value = V[2]
        end
        if n == 0
            thetaval2 = V[1]
            Θ = (thetaval1 - thetaval2) / (4.0 * Δt)
            return (value, Δ, Γ, Θ)
        end
    end
end

function crr_fut(optType, F, K, r_csa_mm_pct, vol_pct, dt_trade, dt_expiry, stepsPerDay)
    # no repo, no dividends, no settlement lags

    dt_option_settle = bd(dt_trade)
    dt_stock_settle = bd(dt_trade)

    ndays = (dt_expiry - dt_trade).value
    T::Float64 = ndays / 365
    N = ceil(Int, stepsPerDay * ndays) + 2
    #Δt = T / N
    Δt = T / (N - 2) # we need to extend the tree by 2 steps to get the right delta, but we don't want to increase the stdev on the terminal slice
    r_csa = log(1.0 + 0.01 * r_csa_mm_pct * ndays / 360) / T
    
    σ = vol_pct / 100.0

    dates = dt_trade:Day(1):dt_expiry

    R = exp(r_csa * Δt)

    U = exp(σ * √Δt)
    D = 1.0 / U
    p = (1.0 - D) / (U - D)
    q = 1.0 - p

    local V::Vector{Float64}
    local px::Float64
    local x::Float64
    local y::Float64
    local dt::Date
    local mult::Float64
    #local F::Float64
    local U::Float64
    local D::Float64
    
    if lowercase(optType)[1] == 'c'
        mult = 1.0
    else
        mult = -1.0
    end

    thetaval1 = thetaval2 = Θ = value = Δ = Γ = 0.0
    V = [max(0.0, mult * (F * U^i * D^(N - i) - K)) for i = 0:N]
    for n = N-1:-1:0
        dt = dt_trade + Day(ceil((n - 2) / (N - 2) * ndays))
        if dt == bd(dt)
            for i = 0:n
                y = (q * V[i+1] + p * V[i+2]) / R
                px = F * U^i * D^(n - i)
                x = mult * (px - K)
                V[i+1] = max(x, y)
            end
        else
            for i = 0:n
                V[i+1] = (q * V[i+1] + p * V[i+2]) / R
            end
        end
        if n == 4
            thetaval1 = V[3]
        end
        if n == 2
            Δ = (V[3] - V[1]) / (F * U^2 - F * D^2)
            Γ = 2 * ((V[3] - V[2]) / (F * U^2 - F * U * D) - (V[2] - V[1]) / (F * U * D - F * D^2)) / (F * U^2 - F * D^2)
            value = V[2]
        end
        if n == 0
            thetaval2 = V[1]
            Θ = (thetaval1 - thetaval2) / (4.0 * Δt)
            return (value, Δ, Γ, Θ)
        end
    end
end

function lr(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)

    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end
    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)

    ndays = (dt_expiry - dt_trade).value
    T = ndays / 365
    N = ceil(Int, stepsPerDay * ndays) + 2
    if iseven(N)
        N::Int = N + 1
    end
    Δt = T / (N - 2)
    r_csa = log(1.0 + 0.01 * r_csa_mm_pct * ndays / 360) / T
    r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / T

    σ = vol_pct / 100.0

    dates = dt_trade:Day(1):dt_expiry

    S0 = S - sum(divamt .* (dt_exdiv .> dt_trade) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt_trade) ./ 360))
    σ = σ * S / S0

    R_csa = exp(r_csa * Δt)
    R_repo = exp(r_repo * Δt)

    #=
    	U = exp(σ * √Δt)
    	D = 1 / U
    	p = (R_repo - D) / (U - D)
    	q = 1-p
    	=#

    #==#

    S_adj = S0 * exp((r_repo - r_csa) * T)
    d1 = (log(S_adj / K) + T * (r_csa + σ^2 / 2)) / (σ * √T)
    d2 = d1 - σ * √T
    function ppif(z, n)
        1 / 2 + sign(z) / 2 * sqrt(1 - exp(-(n + 1 / 6) * (z / (n + 1 / 3 + 0.1 / (n + 1)))^2))
    end
    p = ppif(d2, N - 2)
    q = 1 - p
    p2 = ppif(d1, N - 2)
    U = exp(r_repo * Δt) * p2 / p
    D = exp(r_repo * Δt) * (1 - p2) / q
    #==#
    # use pnow,qnow for dt < dt_stock_settle
    pnow = (1.0 - D) / (U - D)
    qnow = 1.0 - pnow
    # use Rnow (i.e. no discounting) for slices with dt < dt_opt_settle
    Rnow = 1.0

    local Z_csa::Float64
    local Z_repo::Float64
    local V::Vector{Float64}
    local p2use::Float64
    local q2use::Float64
    local R2use::Float64
    local px::Float64
    local x::Float64
    local y::Float64
    local tau_csa::Float64
    local tau_repo::Float64
    local dt::Date
    local mult::Float64

    if lowercase(optType)[1] == 'c'
        mult = 1.0
    else
        mult = -1.0
    end

    thetaval1 = thetaval2 = Θ = value = Δ = Γ = 0
    S0 = S0 / (U * D)

    dt = dt_expiry
    dt_exercise_settle = workday(dt, exercise_settlement)
    tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
    tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
    Z_csa = exp(-r_csa * tau_csa)
    Z_repo = exp(-r_repo * tau_repo)
    V = Z_csa * [max(0.0, mult * (S0 * U^i * D^(N - i) / Z_repo - K)) for i = 0:N]
    for n = N-1:-1:0
        dt = dt_trade + Day(ceil((n - 2) / (N - 2) * ndays))
        dt_next = dt_trade + Day(ceil((n + 1 - 2) / (N - 2) * ndays))
        dt_exercise_settle = workday(dt, exercise_settlement)
        tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
        tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
        Z_csa = exp(-r_csa * tau_csa)
        Z_repo = exp(-r_repo * tau_repo)
        if dt_next <= dt_stock_settle
            p2use = pnow
        else
            p2use = p
        end
        q2use = 1.0 - p2use
        if dt_next <= dt_option_settle
            R2use = Rnow
        else
            R2use = R_csa
        end
        if dt == bd(dt)
            for i = 0:n
                y = (q2use * V[i+1] + p2use * V[i+2]) / R2use
                px = S0 * U^i * D^(n - i) + sum(divamt .* (dt_exdiv .> dt) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt) ./ 360))
                x = Z_csa * mult * (px / Z_repo - K)
                V[i+1] = max(x, y)
            end
        else
            for i = 0:n
                V[i+1] = (q2use * V[i+1] + p2use * V[i+2]) / R2use
            end
        end
        if n == 4
            thetaval1 = V[3]
        end
        if n == 2
            Δ = (V[3] - V[1]) / (S0 * U^2 - S0 * D^2)
            Γ = 2 * ((V[3] - V[2]) / (S0 * U^2 - S0 * U * D) - (V[2] - V[1]) / (S0 * U * D - S0 * D^2)) / (S0 * U^2 - S0 * D^2)
            value = V[2]
        end
        if n == 0
            thetaval2 = V[1]
            Θ = (thetaval1 - thetaval2) / (4 * Δt)
            return (value, Δ, Γ, Θ)
        end
    end
    return V[1]
end

function trinomial(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay, λ=sqrt(2pi) / 2)

    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end
    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)

    ndays = (dt_expiry - dt_trade).value
    T = ndays / 365
    N = ceil(Int, stepsPerDay * ndays) + 1
    Δt = T / (N - 1)
    r_csa = log(1.0 + 0.01 * r_csa_mm_pct * ndays / 360) / T
    r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / T

    σ = vol_pct / 100

    dates = dt_trade:Day(1):dt_expiry
    S0 = S - sum(divamt .* (dt_exdiv .> dt_trade) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt_trade) ./ 360))
    σ = σ * S / S0

    R_csa = exp(r_csa * Δt)
    R_repo = exp(r_repo * Δt)

    ##
    u = exp(λ * σ * sqrt(Δt))
    m = 1
    d = 1 / u
    q = r_csa - r_repo  ## dummy for the div yld - we have discrete divs
    M = exp(r_repo * Δt)
    Σ = M^2 * (-1 + exp(σ^2 * Δt))
    pᵤ = (u * (Σ + M^2 - M) - (M - 1)) / ((u - 1) * (u^2 - 1))
    pdown = (u^2 * (Σ + M^2 - M) - u^3 * (M - 1)) / ((u - 1) * (u^2 - 1))
    pₘ = 1 - pᵤ - pdown
    #pᵤ = 1 / 6 + (r_repo - σ^2 / 2) * sqrt(Δt / (12 * σ^2))
    #pₘ = 2 / 3
    # pdown = 1-pm - pu
    #p_u_now = 1 / 6 + (0 - σ^2 / 2) * sqrt(Δt / (12 * σ^2)) #for dt < dt_stock_settle
    p_u_now = (u * (Σ + 1^2 - 1) - (1 - 1)) / ((u - 1) * (u^2 - 1))
    p_d_now = (u^2 * (Σ + 1^2 - 1) - u^3 * (1 - 1)) / ((u - 1) * (u^2 - 1))
    p_m_now = 1 - p_u_now - p_d_now
    ##
    # use Rnow (i.e. no discounting) for slices with dt < dt_opt_settle
    Rnow = 1.0
    local Z_csa::Float64
    local Z_repo::Float64
    local V::Vector{Float64}
    local p2use::Float64
    local q2use::Float64
    local R2use::Float64
    local px::Float64
    local x::Float64
    local y::Float64
    local tau_csa::Float64
    local tau_repo::Float64
    local dt::Date
    local mult::Float64

    if lowercase(optType)[1] == 'c'
        mult = 1.0
    else
        mult = -1.0
    end

    thetaval1 = thetaval2 = Θ = value = Δ = Γ = 0.0

    dt = dt_expiry
    dt_exercise_settle = workday(dt, exercise_settlement)
    tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
    tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
    Z_csa = exp(-r_csa * tau_csa)
    Z_repo = exp(-r_repo * tau_repo)
    V = Z_csa * [max(0.0, mult * (S0 * u^i / Z_repo - K)) for i = -N:N]
    for n = N-1:-1:0
        dt = dt_trade + Day(ceil((n - 1) / (N - 1) * ndays))
        dt_next = dt_trade + Day(ceil((n + 1 - 1) / (N - 1) * ndays))
        #=
        dt_exercise_settle = workday(dt, exercise_settlement)
        tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
        tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
        Z_csa = exp(-r_csa * tau_csa)
        Z_repo = exp(-r_repo * tau_repo)
        =#
        if dt_next <= dt_stock_settle
            p_u = p_u_now
            p_d = p_d_now
            p_m = p_m_now
        else
            p_u = pᵤ
            p_d = pdown
            p_m = pₘ
        end
        if dt_next <= dt_option_settle
            R = Rnow
        else
            R = R_csa
        end
        if dt == bd(dt)
            for i = -n:n
                y = (p_d * V[i+n+1] + p_m * V[i+n+2] + p_u * V[i+n+3]) / R
                px = S0 * u^i + sum(divamt .* (dt_exdiv .> dt) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt) ./ 360))
                x = Z_csa * mult * (px / Z_repo - K)
                V[i+n+1] = max(x, y)
            end
        else
            for i = -n:n
                V[i+n+1] = (p_d * V[i+n+1] + p_m * V[i+n+2] + p_u * V[i+n+3]) / R
            end
        end
        if n == 2
            thetaval1 = V[3]
        end
        if n == 1
            Δ = (V[3] - V[1]) / (S0 * u - S0 * d)
            Γ = 2 * ((V[3] - V[2]) / (S0 * u - S0 * m) - (V[2] - V[1]) / (S0 * m - S0 * d)) / (S0 * u - S0 * d)
            value = V[2]
        end
        if n == 0
            thetaval2 = V[1]
            Θ = (thetaval1 - thetaval2) / (2 * Δt)
            return (value, Δ, Γ, Θ)
        end
    end
end

function option_am(treeType, optType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    if lowercase(treeType) == "crr"
        fn = crr
    elseif lowercase(treeType) == "lr"
        fn = lr
    elseif lowercase(treeType) == "trinomial"
        fn = trinomial
    end
    return fn(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
end

function implied_vol_am(treeType, optPrice, optType, S, K, r_csa_mm_pct, r_repo_mm_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    if lowercase(treeType) == "crr"
        fn = crr
    elseif lowercase(treeType) == "lr"
        fn = lr
    elseif lowercase(treeType) == "trinomial"
        fn = trinomial
    end
    obj(x) = (fn(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, x, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1] - optPrice)^2
    res = Optim.optimize(obj, 0, 1000)
    vol = Optim.minimizer(res)
    vol
end

function implied_vol_am_fut(treeType, optPrice, optType, F, K, r_csa_mm_pct, dt_trade, dt_expiry, stepsPerDay)
    if lowercase(treeType) == "crr_fut"
        fn = crr_fut
    elseif lowercase(treeType) == "lr_fut"
        fn = lr_fut
    elseif lowercase(treeType) == "trinomial_fut"
        fn = trinomial_fut
    end
    obj(x) = (fn(optType, F, K, r_csa_mm_pct, x, dt_trade, dt_expiry, stepsPerDay)[1] - optPrice)^2
    res = Optim.optimize(obj, 0, 1000)
    vol = Optim.minimizer(res)
    vol
end

function implied_vol_and_repo_from_put_and_call_am(treeType, putPrice, callPrice, S, K, r_csa_mm_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    # this function simultaneously solves for the implied volatility and repo rate for a pair of same-strike options at a given csa rate
    # VERY USEFUL if you don't have a good guess for the repo rate from your Prime Broker
    if lowercase(treeType) == "crr"
        fn = crr
    elseif lowercase(treeType) == "lr"
        fn = lr
    elseif lowercase(treeType) == "trinomial"
        fn = trinomial
    end
    obj(x) = (callPrice - fn("call", S, K, r_csa_mm_pct, x[1], x[2]^2, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1])^2 +
             (putPrice - fn("put", S, K, r_csa_mm_pct, x[1], x[2]^2, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1])^2

    initial_x = [r_csa_mm_pct, sqrt(50)]
    res = Optim.optimize(obj, initial_x)
    if Optim.converged(res)
        r_repo_mm_pct = Optim.minimizer(res)[1]
        vol_pct = Optim.minimizer(res)[2]^2
        return (r_repo_mm_pct, vol_pct)
        #println("call and put are both correctly priced when:")
        #println("r_repo_mm_pct = $(round(r_repo_mm_pct,digits=2)) and")
        #println("vol_pct = $(round(vol_pct))\n")
    else
        return (NaN, NaN)
        #println("Optimizer did not find a solution")
    end
end

function binomial_exercise_simulator(optType, nodes, S0, r_csa, r_repo, σ, divamt, dt_exdiv, dt_divpmt, U, D, p_up, dt_trade, dt_expiry, nSim) #p = p_up, q = p_down = 1-p
    # THIS FUNCTION ASSUMES stepsPerDay = 1
    sort!(nodes, :date)
    dates = dt_trade:Day(1):dt_expiry
    ndays = length(dates)
    divfvs = zeros(length(dates))
    for i in 1:ndays
        divfvs[i] = sum(divamt .* (dt_exdiv .> dates[i]) .* exp.(-r_csa .* Dates.value.(dt_divpmt .- dates[i]) ./ 360))
    end

    if nrow(nodes) == 0
        println("EARLY EXERCISE CONDITION NOT MET ANYWHERE IN TREE")
        return []
    else
        # define a dataframe that summarizes the exercise info from nodes
        # to make code more readable we just call it X
        X = DataFrame(sliceNum=0:ndays-1, date=dates, exPrice=zeros(ndays) .+ NaN, exPriceSigs=zeros(ndays) .+ NaN, sExPrice=zeros(ndays) .+ NaN)
        Threads.@threads for r in 1:nrow(nodes)
            i = findall(dates .== nodes.date[r])[1]
            X.exPrice[i] = nodes.exPrice[r]
        end
        # recall that Sₜ=S₀exp((μ-σ²/2)t + σϵ√t) where ϵ is drawn from N(0,1) implies <Sₜ> = S₀exp(μt)  and var(Sₜ) = S₀^2 * exp(2μt) * (-1 + exp(σ²t))
        # so pₜ is "n sigmas" away from <Sₜ> where n is defined by n = (pₜ - <Sₜ>) / sqrt(var(Sₜ)) = [pₜ/S₀ - exp(μt)] / [exp(μt) * sqrt(-1 + exp(σ²t))]
        Threads.@threads for r = 1:nrow(X)
            if X.exPrice[r] > 0
                t = X.sliceNum[r] / 360
                X.exPriceSigs[r] = abs(((X.exPrice[r] - divfvs[r]) / S0 - exp(r_repo * t)) / (exp(r_repo * t) * sqrt(-1 + exp(σ^2 * t))))
            end
        end
        indx = findall(X.exPriceSigs .>= 10)
        if length(indx) > 0
            X.exPrice[indx] .= NaN
        end
        if sum(divamt) == 0
            indx = findall(X.exPrice .> 0)
            # fit a 4th order polynomial thru exercise nodes
            if length(indx) > 0
                yfit = ls_fit(X.sliceNum[indx], X.exPrice[indx], X.sliceNum[indx[1]]:X.sliceNum[indx[end]], 4)[2]
                X.sExPrice[indx[1]:indx[end]] = yfit
            end
        else
            i_exdiv = zeros(Int, length(dt_exdiv))
            for i = 1:length(dt_exdiv)
                i_exdiv[i] = findall(X.date .== dt_exdiv[i])[1]
            end
            indx = findall(X.exPrice .> 0 .&& X.date .< dt_exdiv[1])
            if length(indx) > 0
                #yfit = ls_fit(X.sliceNum[indx], X.exPrice[indx], X.sliceNum[indx[1]]:i_exdiv[1]-2, 4)[2]
                #X.sExPrice[indx[1]:i_exdiv[1]-1] = yfit
                yfit = ls_fit(X.sliceNum[indx], X.exPrice[indx], X.sliceNum[indx[1]]:X.sliceNum[indx[end]], 4)[2]
                X.sExPrice[indx[1]:indx[end]] = yfit
            end
            for i = 2:length(divamt)
                indx = findall(X.exPrice .> 0 .&& X.date .>= dt_exdiv[i-1] .&& X.date .< dt_exdiv[i])
                if length(indx) > 0
                    #yfit = ls_fit(X.sliceNum[indx], X.exPrice[indx], X.sliceNum[indx[1]]:i_exdiv[i]-2, 4)[2]
                    #X.sExPrice[indx[1]:i_exdiv[i]-1] = yfit
                    yfit = ls_fit(X.sliceNum[indx], X.exPrice[indx], X.sliceNum[indx[1]]:X.sliceNum[indx[end]], 4)[2]
                    X.sExPrice[indx[1]:indx[end]] = yfit
                end
            end
            indx = findall(X.exPrice .> 0 .&& X.date .>= dt_exdiv[end])
            if length(indx) > 0
                #yfit = ls_fit(X.sliceNum[indx], X.exPrice[indx], X.sliceNum[indx[1]]:X.sliceNum[end], 4)[2]
                #X.sExPrice[indx[1]:end] = yfit
                yfit = ls_fit(X.sliceNum[indx], X.exPrice[indx], X.sliceNum[indx[1]]:X.sliceNum[indx[end]], 4)[2]
                X.sExPrice[indx[1]:indx[end]] = yfit
            end
        end


        #=

        let's map the smoothed exercise boundary value at slice number n to a number of up moves, i.e. 
        S0*U^i*D^(n-iM) + divfvs[n] = yhat[n] at slice n
        U^i D^(n-i) = (yhat[n]-divfvs[n])/S0
        i*logU + (n-i)logD = log((yhat[n]-divfvs[n])/S0)
        i(logU-logD) = -nlogD + log((yhat[n]-divfvs[n])/S0)
        i = log((yhat[n]-divfvs[n])/S0)-nlogD)/(logU-logD) #any inequality sign stays the same because logU > logD

        to compute probability of crossing the smoothed call boundary, we need
        i >= iMin = ceil(Int,log((yhat[n]-divfvs[n])/S0)-nlogD)/(logU-logD) #inequality sign stays >= because logU > logD

        in other words, if on a given path that hasn't had early exercise before slice n, a call will be exercised if the cumulative number of up 
        moves by slice n is >= iMin

        conversely, to compute probability of crossing the smoothed put boundary, we need
        i <= iMax = floor(Int,log((yhat[n]-divfvs[n])/S0)-nlogD)/(logU-logD)

        in other words, if on a given path that hasn't had early exercise before slice n, a put will be exercised if the cumulative number of up 
        moves by slice n is <= iMax
        =#

        d = Binomial(1, p_up)
        # ndays rows, nSim columns, first row should be filled with 0s
        sims = rand(d, ndays, Int(nSim)) # here 0s indicate downs and 1s indicate ups
        sims[1, :] .= 0
        numUps = cumsum(sims, dims=1)

        if lowercase(optType)[1] == 'c'
            v = copy(X.sExPrice)
            v[isnan.(v)] .= S0 * U^(ndays + 10) + divfvs[1]
            iMin = ceil.(Int, (log.(max.(0.0001, v - divfvs) ./ S0) - X.sliceNum .* log(D)) ./ (log(U) - log(D)))
            exCond = (numUps .>= iMin)
        elseif lowercase(optType)[1] == 'p'
            v = copy(X.sExPrice)
            v[isnan.(v)] .= 0.0001
            iMax = floor.(Int, (log.(max.(0.0001, v - divfvs) ./ S0) - X.sliceNum .* log(D)) ./ (log(U) - log(D)))
            exCond = (numUps .<= iMax)
        end
        Threads.@threads for c = 1:size(exCond, 2)
            local indx = findall(exCond[:, c])
            if length(indx) > 0
                exCond[indx[1]:end, c] .= true
            end
        end
        cumExProb = zeros(ndays)
        Threads.@threads for i = 1:ndays
            cumExProb[i] = length(findall(exCond[i, :])) / nSim
        end
        X.cumExProb = cumExProb
        return X
    end
end

function crr_extra_output(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay=1)

    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end

    boundaryThresh = 0

    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)

    ndays = (dt_expiry - dt_trade).value
    T = ndays / 365
    N = ceil(Int, stepsPerDay * ndays)
    Δt = T / N
    r_csa = log(1.0 + 0.01 * r_csa_mm_pct * ndays / 360) / T
    r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / T

    σ = vol_pct / 100

    dates = dt_trade:Day(1):dt_expiry


    # we will discount using R_csa from i_option_settle to end
    # we will grow S0 an R_repo from i_stock_settle to end
    # this means we'll go to the LAST slice on dt_option_settle and dt_stock_settle
    #i_option_settle = 1 + stepsPerDay * (dt_option_settle - dt_trade).value
    #i_stock_settle = 1 + stepsPerDay * (dt_stock_settle - dt_trade).value
    S0 = S - sum(divamt .* (dt_exdiv .> dt_trade) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt_trade) ./ 360))
    σ = σ * S / S0

    R_csa = exp(r_csa * Δt)
    R_repo = exp(r_repo * Δt)

    U = exp(σ * √Δt)
    D = 1.0 / U
    p = (R_repo - D) / (U - D)
    q = 1.0 - p
    # use pnow,qnow for dt < dt_stock_settle
    pnow = (1.0 - D) / (U - D)
    qnow = 1.0 - pnow
    # use Rnow (i.e. no discounting) for slices with dt < dt_opt_settle
    Rnow = 1.0

    local Z_csa::Float64
    local Z_repo::Float64
    local V::Vector{Float64}
    local p2use::Float64
    local q2use::Float64
    local R2use::Float64
    local px::Float64
    local x::Float64
    local y::Float64
    local tau_csa::Float64
    local tau_repo::Float64
    local dt::Date

    nodes = DataFrame(date=Date[], exPrice=Float64[])
    if lowercase(optType)[1] == 'p'
        dt = dt_expiry
        dt_exercise_settle = workday(dt, exercise_settlement)
        tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
        tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
        Z_csa = exp(-r_csa * tau_csa)
        Z_repo = exp(-r_repo * tau_repo)
        V = Z_csa * [max(0.0, K - S0 * U^i * D^(N - i) / Z_repo) for i = 0:N]
        for n = N-1:-1:0
            dt = dt_trade + Day(ceil(n / N * ndays))
            dt_next = dt_trade + Day(ceil((n + 1) / N * ndays))
            dt_exercise_settle = workday(dt, exercise_settlement)
            tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
            tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
            Z_csa = exp(-r_csa * tau_csa)
            Z_repo = exp(-r_repo * tau_repo)
            if dt_next <= dt_stock_settle
                p2use = pnow
            else
                p2use = p
            end
            q2use = 1.0 - p2use
            if dt_next <= dt_option_settle
                R2use = Rnow
            else
                R2use = R_csa
            end
            if dt == bd(dt)
                for i = 0:n
                    y = (q2use * V[i+1] + p2use * V[i+2]) / R2use
                    px = S0 * U^i * D^(n - i) + sum(divamt .* (dt_exdiv .> dt) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt) ./ 360))
                    x = Z_csa * (K - px / Z_repo)
                    V[i+1] = max(x, y)
                    if x > y + boundaryThresh
                        if nrow(nodes) == 0 || dt < nodes.date[end]
                            push!(nodes, [dt, px])
                        elseif nrow(nodes) > 0 && dt == nodes.date[end] #keep highest put early ex px if multiple slices per day
                            # this logic is put-specific since i=0:n is in strictly
                            # increasing price order and we want the highest early
                            # ex price for the put boundary
                            if px > nodes.exPrice[end]
                                pop!(nodes)
                                push!(nodes, [dt, px])
                            end
                        end
                    end
                end
            else
                for i = 0:n
                    V[i+1] = (q2use * V[i+1] + p2use * V[i+2]) / R2use
                end
            end
        end
    elseif lowercase(optType)[1] == 'c'
        dt = dt_expiry
        dt_exercise_settle = workday(dt, exercise_settlement)
        tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
        tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
        Z_csa = exp(-r_csa * tau_csa)
        Z_repo = exp(-r_repo * tau_repo)
        V = Z_csa * [max(0.0, S0 * U^i * D^(N - i) / Z_repo - K) for i = 0:N]
        for n = N-1:-1:0
            dt = dt_trade + Day(ceil(n / N * ndays))
            dt_next = dt_trade + Day(ceil((n + 1) / N * ndays))
            dt_exercise_settle = workday(dt, exercise_settlement)
            tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
            tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
            Z_csa = exp(-r_csa * tau_csa)
            Z_repo = exp(-r_repo * tau_repo)
            if dt_next <= dt_stock_settle
                p2use = pnow
            else
                p2use = p
            end
            q2use = 1.0 - p2use
            if dt_next <= dt_option_settle
                R2use = Rnow
            else
                R2use = R_csa
            end
            if dt == bd(dt)
                for i = 0:n
                    y = (q2use * V[i+1] + p2use * V[i+2]) / R2use
                    px = S0 * U^i * D^(n - i) + sum(divamt .* (dt_exdiv .> dt) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt) ./ 360))
                    x = Z_csa * (px / Z_repo - K)
                    V[i+1] = max(x, y)
                    # this logic is call-specific since i=0:n is in strictly
                    # increasing price order and we want the lowest early
                    # ex price for the call boundary will be the first encountered
                    # as i runs from 0 to n
                    if x > y + boundaryThresh
                        if nrow(nodes) == 0 || dt < nodes.date[end]
                            push!(nodes, [dt, px])
                        elseif nrow(nodes) > 0 && dt == nodes.date[end] #want lowest call ex px on day if stepsPerDay > 1
                            if px < nodes.exPrice[end]
                                pop!(nodes)
                                push!(nodes, [dt, px])
                            end
                        end
                    end
                end
            else
                for i = 0:n
                    V[i+1] = (q2use * V[i+1] + p2use * V[i+2]) / R2use
                end
            end
        end
    end
    X = binomial_exercise_simulator(optType, nodes, S0, r_csa, r_repo, σ, divamt, dt_exdiv, dt_divpmt, U, D, p, dt_trade, dt_expiry, 10000)
    return (V[1], X, optType)
    #return (V[1], nodes, optType)
end

function lr_extra_output(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay=1)

    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end

    boundaryThresh = 0

    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)

    ndays = (dt_expiry - dt_trade).value
    T = ndays / 365
    N = ceil(Int, stepsPerDay * ndays)
    if N % 2 < 0.5
        N::Int = N + 1
    end
    Δt = T / N
    r_csa = log(1.0 + 0.01 * r_csa_mm_pct * ndays / 360) / T
    r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / T

    σ = vol_pct / 100.0

    dates = dt_trade:Day(1):dt_expiry

    S0 = S - sum(divamt .* (dt_exdiv .> dt_trade) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt_trade) ./ 360))
    σ = σ * S / S0

    R_csa = exp(r_csa * Δt)
    R_repo = exp(r_repo * Δt)

    #=
    	U = exp(σ * √Δt)
    	D = 1 / U
    	p = (R_repo - D) / (U - D)
    	q = 1-p
    	=#

    #==#

    S_adj = S0 * exp((r_repo - r_csa) * T)
    d1 = (log(S_adj / K) + T * (r_csa + σ^2 / 2)) / (σ * √T)
    d2 = d1 - σ * √T
    function ppif(z, n)
        1 / 2 + sign(z) / 2 * sqrt(1 - exp(-(n + 1 / 6) * (z / (n + 1 / 3 + 0.1 / (n + 1)))^2))
    end
    p = ppif(d2, N)
    q = 1 - p
    p2 = ppif(d1, N)
    U = exp(r_repo * Δt) * p2 / p
    D = exp(r_repo * Δt) * (1 - p2) / q
    #==#
    # use pnow,qnow for dt < dt_stock_settle
    pnow = (1.0 - D) / (U - D)
    qnow = 1.0 - pnow
    # use Rnow (i.e. no discounting) for slices with dt < dt_opt_settle
    Rnow = 1.0

    local Z_csa::Float64
    local Z_repo::Float64
    local V::Vector{Float64}
    local p2use::Float64
    local q2use::Float64
    local R2use::Float64
    local px::Float64
    local x::Float64
    local y::Float64
    local tau_csa::Float64
    local tau_repo::Float64
    local dt::Date

    nodes = DataFrame(date=Date[], exPrice=Float64[])
    if lowercase(optType)[1] == 'p'
        dt = dt_expiry
        dt_exercise_settle = workday(dt, exercise_settlement)
        tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
        tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
        Z_csa = exp(-r_csa * tau_csa)
        Z_repo = exp(-r_repo * tau_repo)
        V = Z_csa * [max(0.0, K - S0 * U^i * D^(N - i) / Z_repo) for i = 0:N]
        for n = N-1:-1:0
            dt = dt_trade + Day(ceil(n / N * ndays))
            dt_next = dt_trade + Day(ceil((n + 1) / N * ndays))
            dt_exercise_settle = workday(dt, exercise_settlement)
            tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
            tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
            Z_csa = exp(-r_csa * tau_csa)
            Z_repo = exp(-r_repo * tau_repo)
            if dt_next <= dt_stock_settle
                p2use = pnow
            else
                p2use = p
            end
            q2use = 1.0 - p2use
            if dt_next <= dt_option_settle
                R2use = Rnow
            else
                R2use = R_csa
            end
            if dt == bd(dt)
                for i = 0:n
                    y = (q2use * V[i+1] + p2use * V[i+2]) / R2use
                    px = S0 * U^i * D^(n - i) + sum(divamt .* (dt_exdiv .> dt) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt) ./ 360))
                    x = Z_csa * (K - px / Z_repo)
                    V[i+1] = max(x, y)
                    if x > y + boundaryThresh
                        if nrow(nodes) == 0 || dt < nodes.date[end]
                            push!(nodes, [dt, px])
                        elseif nrow(nodes) > 0 && dt == nodes.date[end] #keep highest put early ex px if multiple slices per day
                            # this logic is put-specific since i=0:n is in strictly
                            # increasing price order and we want the highest early
                            # ex price for the put boundary
                            if px > nodes.exPrice[end]
                                pop!(nodes)
                                push!(nodes, [dt, px])
                            end
                        end
                    end
                end
            else
                for i = 0:n
                    V[i+1] = (q2use * V[i+1] + p2use * V[i+2]) / R2use
                end
            end
        end
    elseif lowercase(optType)[1] == 'c'
        dt = dt_expiry
        dt_exercise_settle = workday(dt, exercise_settlement)
        tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
        tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
        Z_csa = exp(-r_csa * tau_csa)
        Z_repo = exp(-r_repo * tau_repo)
        V = Z_csa * [max(0.0, S0 * U^i * D^(N - i) / Z_repo - K) for i = 0:N]

        for n = N-1:-1:0
            dt = dt_trade + Day(ceil(n / N * ndays))
            dt_next = dt_trade + Day(ceil((n + 1) / N * ndays))
            dt_exercise_settle = workday(dt, exercise_settlement)
            tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
            tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
            Z_csa = exp(-r_csa * tau_csa)
            Z_repo = exp(-r_repo * tau_repo)
            if dt_next <= dt_stock_settle
                p2use = pnow
            else
                p2use = p
            end
            q2use = 1.0 - p2use
            if dt_next <= dt_option_settle
                R2use = Rnow
            else
                R2use = R_csa
            end
            if dt == bd(dt)
                for i = 0:n
                    y = (q2use * V[i+1] + p2use * V[i+2]) / R2use
                    px = S0 * U^i * D^(n - i) + sum(divamt .* (dt_exdiv .> dt) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt) ./ 360))
                    x = Z_csa * (px / Z_repo - K)
                    V[i+1] = max(x, y)
                    # this logic is call-specific since i=0:n is in strictly
                    # increasing price order and we want the lowest early
                    # ex price for the call boundary will be the first encountered
                    # as i runs from 0 to n
                    if x > y + boundaryThresh
                        if nrow(nodes) == 0 || dt < nodes.date[end]
                            push!(nodes, [dt, px])
                        elseif nrow(nodes) > 0 && dt == nodes.date[end] #want lowest call ex px on day if stepsPerDay > 1
                            if px < nodes.exPrice[end]
                                pop!(nodes)
                                push!(nodes, [dt, px])
                            end
                        end
                    end
                end
            else
                for i = 0:n
                    V[i+1] = (q2use * V[i+1] + p2use * V[i+2]) / R2use
                end
            end
        end
    end
    X = binomial_exercise_simulator(optType, nodes, S0, r_csa, r_repo, σ, divamt, dt_exdiv, dt_divpmt, U, D, p, dt_trade, dt_expiry, 10000)
    return (V[1], X, optType)
end

function trinomial_exercise_simulator(optType, nodes, S0, r_csa, r_repo, σ, divamt, dt_exdiv, dt_divpmt, u, p_u, dt_trade, dt_expiry, nSim) #p = p_up, q = p_down = 1-p
    # THIS FUNCTION ASSUMES stepsPerDay = 1
    d = 1 / u
    m = 1
    p_m = 2 / 3
    p_d = 1 - p_u - p_m

    sort!(nodes, :date)
    dates = dt_trade:Day(1):dt_expiry
    ndays = length(dates)
    divfvs = zeros(length(dates))
    for i in 1:ndays
        divfvs[i] = sum(divamt .* (dt_exdiv .> dates[i]) .* exp.(-r_csa .* Dates.value.(dt_divpmt .- dates[i]) ./ 360))
    end

    if nrow(nodes) == 0
        println("EARLY EXERCISE CONDITION NOT MET ANYWHERE IN TREE")
        return []
    else
        # define a dataframe that summarizes the exercise info from nodes
        # to make code more readable we just call it X
        X = DataFrame(sliceNum=0:ndays-1, date=dates, exPrice=zeros(ndays) .+ NaN, exPriceSigs=zeros(ndays) .+ NaN, sExPrice=zeros(ndays) .+ NaN)
        Threads.@threads for r in 1:nrow(nodes)
            i = findall(dates .== nodes.date[r])[1]
            X.exPrice[i] = nodes.exPrice[r]
        end
        # recall that Sₜ=S₀exp((μ-σ²/2)t + σϵ√t) where ϵ is drawn from N(0,1) implies <Sₜ> = S₀exp(μt)  and var(Sₜ) = S₀^2 * exp(2μt) * (-1 + exp(σ²t))
        # so pₜ is "n sigmas" away from <Sₜ> where n is defined by n = (pₜ - <Sₜ>) / sqrt(var(Sₜ)) = [pₜ/S₀ - exp(μt)] / [exp(μt) * sqrt(-1 + exp(σ²t))]
        Threads.@threads for r = 1:nrow(X)
            if X.exPrice[r] > 0
                t = X.sliceNum[r] / 360
                X.exPriceSigs[r] = abs(((X.exPrice[r] - divfvs[r]) / S0 - exp(r_repo * t)) / (exp(r_repo * t) * sqrt(-1 + exp(σ^2 * t))))
            end
        end
        indx = findall(X.exPriceSigs .>= 10)
        if length(indx) > 0
            X.exPrice[indx] .= NaN
        end
        if sum(divamt) == 0
            indx = findall(X.exPrice .> 0)
            # fit a 4th order polynomial thru exercise nodes
            yfit = ls_fit(X.sliceNum[indx], X.exPrice[indx], X.sliceNum[indx[1]]:X.sliceNum[indx[end]], 4)[2]
            X.sExPrice[indx[1]:indx[end]] = yfit
        else
            indx = findall(X.exPrice .> 0 .&& X.date .< dt_exdiv[1])
            if length(indx) > 0
                yfit = ls_fit(X.sliceNum[indx], X.exPrice[indx], X.sliceNum[indx[1]]:X.sliceNum[indx[end]], 4)[2]
                X.sExPrice[indx[1]:indx[end]] = yfit
            end
            for i = 2:length(divamt)
                indx = findall(X.exPrice .> 0 .&& X.date .>= dt_exdiv[i-1] .&& X.date .< dt_exdiv[i])
                if length(indx) > 0
                    yfit = ls_fit(X.sliceNum[indx], X.exPrice[indx], X.sliceNum[indx[1]]:X.sliceNum[indx[end]], 4)[2]
                    X.sExPrice[indx[1]:indx[end]] = yfit
                end
            end
            indx = findall(X.exPrice .> 0 .&& X.date .>= dt_exdiv[end])
            if length(indx) > 0
                yfit = ls_fit(X.sliceNum[indx], X.exPrice[indx], X.sliceNum[indx[1]]:X.sliceNum[indx[end]], 4)[2]
                X.sExPrice[indx[1]:indx[end]] = yfit
            end
        end

        #=

        let's map the smoothed exercise boundary value at slice number n to a number of up moves, i.e. 
        S0*u^i + divfvs[n] = yhat[n] at slice n (i can run from -n to n)
        u^i = (yhat[n]-divfvs[n])/S0
        i*logu = log((yhat[n]-divfvs[n])/S0)
        i = log((yhat[n]-divfvs[n])/S0)/logu #any inequality sign stays the same because logu > 0

        to compute probability of crossing the smoothed call boundary, we need
        i >= iMin = ceil(Int,log((yhat[n]-divfvs[n])/S0))/logu #inequality sign stays >= because logu > 0

        in other words, if on a given path that hasn't had early exercise before slice n, a call will be exercised if the cumulative number of up 
        moves by slice n is >= iMin

        conversely, to compute probability of crossing the smoothed put boundary, we need
        i <= iMax = floor(Int,log((yhat[n]-divfvs[n])/S0))/logu

        in other words, if on a given path that hasn't had early exercise before slice n, a put will be exercised if the cumulative number of up 
        moves by slice n is <= iMax
        =#

        d = Multinomial(1, [p_u, 2 / 3, 1 / 3 - p_u])
        # ndays rows, nSim columns, first row should be filled with 0s
        sims = rand(d, ndays, Int(nSim)) # here 0s indicate downs and 1s indicate ups
        Threads.@threads for c = 1:nSim
            sims[1, c] = [0, 1, 0]
        end
        numUps = zeros(ndays, nSim)
        Threads.@threads for c = 1:nSim
            for r = 2:ndays
                numUps[r, c] = numUps[r-1, c]
                if sims[r, c] == [1, 0, 0]
                    numUps[r, c] = numUps[r, c] + 1
                elseif sims[r, c] == [0, 0, 1]
                    numUps[r, c] = numUps[r, c] - 1
                end
            end
        end

        if lowercase(optType)[1] == 'c'
            v = copy(X.sExPrice)
            v[isnan.(v)] .= S0 * u^(ndays + 10) + divfvs[1]
            iMin = ceil.(Int, log.(max.(0.0001, v - divfvs) / S0) / log(u))
            exCond = (numUps .>= iMin)
        elseif lowercase(optType)[1] == 'p'
            v = copy(X.sExPrice)
            v[isnan.(v)] .= 0.0001
            iMax = floor.(Int, log.(max.(0.0001, v - divfvs) / S0) / log(u))
            exCond = (numUps .<= iMax)
        end
        Threads.@threads for c = 1:size(exCond, 2)
            local indx = findall(exCond[:, c])
            if length(indx) > 0
                exCond[indx[1]:end, c] .= true
            end
        end
        cumExProb = zeros(ndays)
        Threads.@threads for i = 1:ndays
            cumExProb[i] = length(findall(exCond[i, :])) / nSim
        end
        X.cumExProb = cumExProb
        return X
    end
end

function trinomial_extra_output(optType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay=1)

    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end

    boundaryThresh = 0

    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)

    ndays = (dt_expiry - dt_trade).value
    T = ndays / 365
    N = ceil(Int, stepsPerDay * ndays)
    Δt = T / N
    r_csa = log(1.0 + 0.01 * r_csa_mm_pct * ndays / 360) / T
    r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / T

    σ = vol_pct / 100

    dates = dt_trade:Day(1):dt_expiry


    # we will discount using R_csa from i_option_settle to end
    # we will grow S0 an R_repo from i_stock_settle to end
    # this means we'll go to the LAST slice on dt_option_settle and dt_stock_settle
    #i_option_settle = 1 + stepsPerDay * (dt_option_settle - dt_trade).value
    #i_stock_settle = 1 + stepsPerDay * (dt_stock_settle - dt_trade).value
    S0 = S - sum(divamt .* (dt_exdiv .> dt_trade) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt_trade) ./ 360))
    σ = σ * S / S0

    R_csa = exp(r_csa * Δt)
    R_repo = exp(r_repo * Δt)

    ##
    u = exp(σ * sqrt(3 * Δt))
    m = 1
    d = 1 / u
    pᵤ = 1 / 6 + (r_repo - σ^2 / 2) * sqrt(Δt / (12 * σ^2))
    pₘ = 2 / 3
    # pdown = 1-pm - pu
    p_u_now = 1 / 6 + (0 - σ^2 / 2) * sqrt(Δt / (12 * σ^2)) #for dt < dt_stock_settle
    ##
    # use Rnow (i.e. no discounting) for slices with dt < dt_opt_settle
    Rnow = 1.0
    local Z_csa::Float64
    local Z_repo::Float64
    local V::Vector{Float64}
    local p2use::Float64
    local q2use::Float64
    local R2use::Float64
    local px::Float64
    local x::Float64
    local y::Float64
    local tau_csa::Float64
    local tau_repo::Float64
    local dt::Date

    nodes = DataFrame(date=Date[], exPrice=Float64[])
    if lowercase(optType)[1] == 'p'
        dt = dt_expiry
        dt_exercise_settle = workday(dt, exercise_settlement)
        tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
        tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
        Z_csa = exp(-r_csa * tau_csa)
        Z_repo = exp(-r_repo * tau_repo)
        V = Z_csa * [max(0.0, K - S0 * u^i / Z_repo) for i = -N:N]
        for n = N-1:-1:0
            dt = dt_trade + Day(ceil(n / N * ndays))
            dt_next = dt_trade + Day(ceil((n + 1) / N * ndays))
            dt_exercise_settle = workday(dt, exercise_settlement)
            tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
            tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
            Z_csa = exp(-r_csa * tau_csa)
            Z_repo = exp(-r_repo * tau_repo)
            if dt_next <= dt_stock_settle
                p_u = p_u_now
            else
                p_u = pᵤ
            end
            if dt_next <= dt_option_settle
                R = Rnow
            else
                R = R_csa
            end
            if dt == bd(dt)
                for i = -n:n
                    y = ((1 - p_u - 2 / 3) * V[i+n+1] + 2 / 3 * V[i+n+2] + p_u * V[i+n+3]) / R
                    px = S0 * u^i + sum(divamt .* (dt_exdiv .> dt) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt) ./ 360))
                    x = Z_csa * (K - px / Z_repo)
                    V[i+n+1] = max(x, y)
                    if x > y + boundaryThresh
                        if nrow(nodes) == 0 || dt < nodes.date[end]
                            push!(nodes, [dt, px])
                        elseif nrow(nodes) > 0 && dt == nodes.date[end] #keep highest put early ex px if multiple slices per day
                            # this logic is put-specific since i=0:n is in strictly
                            # increasing price order and we want the highest early
                            # ex price for the put boundary
                            if px > nodes.exPrice[end]
                                pop!(nodes)
                                push!(nodes, [dt, px])
                            end
                        end
                    end
                end
            else
                for i = -n:n
                    V[i+n+1] = ((1 - p_u - 2 / 3) * V[i+n+1] + 2 / 3 * V[i+n+2] + p_u * V[i+n+3]) / R
                end
            end
        end
    elseif lowercase(optType)[1] == 'c'
        dt = dt_expiry
        dt_exercise_settle = workday(dt, exercise_settlement)
        tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
        tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
        Z_csa = exp(-r_csa * tau_csa)
        Z_repo = exp(-r_repo * tau_repo)
        V = Z_csa * [max(0.0, S0 * u^i / Z_repo - K) for i = -N:N]

        for n = N-1:-1:0
            dt = dt_trade + Day(ceil(n / N * ndays))
            dt_next = dt_trade + Day(ceil((n + 1) / N * ndays))
            dt_exercise_settle = workday(dt, exercise_settlement)
            tau_repo = (dt_exercise_settle - max(dt, dt_stock_settle)).value / 365
            tau_csa = (dt_exercise_settle - max(dt, dt_option_settle)).value / 365
            Z_csa = exp(-r_csa * tau_csa)
            Z_repo = exp(-r_repo * tau_repo)
            if dt_next <= dt_stock_settle
                p_u = p_u_now
            else
                p_u = pᵤ
            end
            if dt_next <= dt_option_settle
                R = Rnow
            else
                R = R_csa
            end
            if dt == bd(dt)
                for i = -n:n
                    y = ((1 - p_u - 2 / 3) * V[i+n+1] + 2 / 3 * V[i+n+2] + p_u * V[i+n+3]) / R
                    px = S0 * u^i + sum(divamt .* (dt_exdiv .> dt) .* exp.(-r_csa * Dates.value.(dt_divpmt - dt) ./ 360))
                    x = Z_csa * (px / Z_repo - K)
                    V[i+n+1] = max(x, y)
                    # this logic is call-specific since i=0:n is in strictly
                    # increasing price order and we want the lowest early
                    # ex price for the call boundary will be the first encountered
                    # as i runs from 0 to n
                    if x > y + boundaryThresh
                        if nrow(nodes) == 0 || dt < nodes.date[end]
                            push!(nodes, [dt, px])
                        elseif nrow(nodes) > 0 && dt == nodes.date[end] #want lowest call ex px on day if stepsPerDay > 1
                            if px < nodes.exPrice[end]
                                pop!(nodes)
                                push!(nodes, [dt, px])
                            end
                        end
                    end
                end
            else
                for i = -n:n
                    V[i+n+1] = ((1 - p_u - 2 / 3) * V[i+n+1] + 2 / 3 * V[i+n+2] + p_u * V[i+n+3]) / R
                end
            end
        end
    end
    X = trinomial_exercise_simulator(optType, nodes, S0, r_csa, r_repo, σ, divamt, dt_exdiv, dt_divpmt, u, pᵤ, dt_trade, dt_expiry, 10000)
    return (V[1], X, optType)
end

function conversion_am(treeType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    if lowercase(treeType) == "crr"
        fn = crr
    elseif lowercase(treeType) == "lr"
        fn = lr
    elseif lowercase(treeType) == "trinomial"
        fn = trinomial
    end

    callPrice = fn("c", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1]
    putPrice = fn("p", S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1]
    conversion = putPrice - callPrice + S - K
    return (conversion)
end

function implied_repo_from_conversion_am(convPrice, treeType, S, K, r_csa_mm_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    #=
     the put/call valuations use r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / T
     so we need to require r_repo_mm_pct > 1 + -36000/ndays
     =#

    T = (dt_expiry - dt_trade).value
    min_repo_pct = 1 - 36000 / T

    global max_allowable_repo_pct = r_csa_mm_pct
    global min_allowable_repo_pct = r_csa_mm_pct

    offset_vec = vcat(30, 0, -1:-1:-3, -10:-10:-100, -200:-100:-900, -1000:-1000:-3000)
    offset_vec = offset_vec[offset_vec.>=(min_repo_pct-r_csa_mm_pct)]
    pvec = zeros(length(offset_vec))
    rvec = r_csa_mm_pct .+ offset_vec
    for i = 1:length(rvec)
        pvec[i] = conversion_am(treeType, S, K, r_csa_mm_pct, rvec[i], vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
        if (i > 1) && ((pvec[i] - convPrice) * (pvec[i-1] - convPrice) < 0) # i.e. the soln lies in between rvec[i] and rvec[i-1]
            max_allowable_repo_pct = rvec[i-1]
            min_allowable_repo_pct = rvec[i]
            break
        elseif (pvec[i] > convPrice)
            return NaN # pvec is a decreasing fn of r so as we walk r in the negative direction, pvec > mktPrice means we can't cross from below and should exit
        end
    end
    if abs(pvec[end]) > 0
        return NaN # pvec[end] will be non-zero only if no soln was found
    end
    obj(x) = (conversion_am(treeType, S, K, r_csa_mm_pct, x, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay) - convPrice)^2
    res = Optim.optimize(obj, min_allowable_repo_pct, max_allowable_repo_pct)
    if Optim.converged(res)
        return Optim.minimizer(res)
    end
end

function implied_dividend_from_conversion_am(convPrice, treeType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, dt_exdiv, dt_divpmt, stepsPerDay, maxDiv2Px)
    n = length(dt_exdiv)
    obj(x) = (conversion_am(treeType, S, K, r_csa_mm_pct, r_repo_mm_pct, vol_pct, dt_trade, dt_expiry, x * ones(n), dt_exdiv[1:n], dt_divpmt[1:n], stepsPerDay) - convPrice)^2
    ans = NaN
    maxrng = 0.01
    while (maxrng < maxDiv2Px)
        res = Optim.optimize(obj, S * (maxrng - 0.01), S * maxrng, abs_tol=0.001)
        if Optim.converged(res) && Optim.minimum(res) < 0.0001
            ans = Optim.minimizer(res)
            break
        else
            maxrng = maxrng + 0.01
        end
    end
    return ans
end

function implied_repo_from_jelly_roll(mktPrice, treeType, S, r_csa_mm_pct, dt_trade, divamt, dt_exdiv, dt_divpmt, stepsPerDay, K1, vol_pct1, dt_expiry1, K2, vol_pct2, dt_expiry2, maxRateInSearchPct=10)
    T1 = (dt_expiry1 - dt_trade).value
    T2 = (dt_expiry2 - dt_trade).value
    min_repo_pct = 1 - 36000 / max(T1, T2)
    #min_repo_pct2 = 1 - 36000 / T2
    #min_repo_pct = max(min_repo_pct1, min_repo_pct2)

    offset_vec =
        vcat(
            maxRateInSearchPct - r_csa_mm_pct,
            -1:-1:-3,
            -10:-10:-100,
            -200:-100:-900, -1000:-1000:-3000
        )
    offset_vec = offset_vec[offset_vec.>=(min_repo_pct-r_csa_mm_pct)]
    pvec = zeros(1 + length(offset_vec))
    rvec = sort(unique(r_csa_mm_pct .+ offset_vec), rev=true)
    for i = 1:length(rvec)
        pvec[i] =
            conversion_am(treeType, S, K2, r_csa_mm_pct, rvec[i], vol_pct2, dt_trade, dt_expiry2, divamt, dt_exdiv, dt_divpmt, stepsPerDay) -
            conversion_am(treeType, S, K1, r_csa_mm_pct, rvec[i], vol_pct1, dt_trade, dt_expiry1, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    end
    global max_allowable_repo = []
    global min_allowable_repo = []
    for i = 2:length(rvec)
        if (pvec[i] - mktPrice) * (pvec[i-1] - mktPrice) < 0 # i.e. the soln lies in between rvec[i] and rvec[i-1]
            global max_allowable_repo = [max_allowable_repo; rvec[i-1]]
            global min_allowable_repo = [min_allowable_repo; rvec[i]]
        end
    end
    global possAns = []
    obj(x) = (conversion_am(treeType, S, K2, r_csa_mm_pct, x, vol_pct2, dt_trade, dt_expiry2, divamt, dt_exdiv, dt_divpmt, stepsPerDay) -
              conversion_am(treeType, S, K1, r_csa_mm_pct, x, vol_pct1, dt_trade, dt_expiry1, divamt, dt_exdiv, dt_divpmt, stepsPerDay) -
              mktPrice)^2
    if length(max_allowable_repo) > 0
        for i = 1:length(max_allowable_repo)
            res = Optim.optimize(
                obj,
                min_allowable_repo[i],
                max_allowable_repo[i],
                abs_tol=1e-5,
            )
            if Optim.converged(res)
                global possAns = [possAns; Optim.minimizer(res)]
            end
        end
    end
    if length(possAns) == 0
        return NaN
    else
        return possAns #[findmin(identity, abs.(possAns .- r_csa))[2]]
    end
end

################ Option pricing functions with yield curve as input

function BS_yc(yc, optType, S, K, r_csa_sprd_pct, r_repo_sprd_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)

    (ycinfo, futs, intervals) = yc
    # note that divamt, dt_exdiv and dt_divpmt need to be array even if there's only one element
    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end
    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)
    dt_exp_settle = workday(dt_expiry, exercise_settlement)

    # put in error trap code in case ycinfo.date[end] < dt_exp_settle

    days = (dt_expiry - dt_trade).value
    t = days / 365
    #r_csa = log(1 + .01 * r_csa_mm_pct* days/360) / t
    #r_repo = log(1 + .01 * r_repo_mm_pct* days/360) / t
    σ = vol_pct / 100
    Z_csa = 1 / (1 + 0.01 * (fwd_rate(dt_option_settle, dt_exp_settle, yc) + r_csa_sprd_pct) * (dt_exp_settle - dt_option_settle).value / 360)

    divPV = 0
    for i = 1:length(divamt)
        divPV = divPV + divamt[i] * df(dt_divpmt[i], yc)
    end
    S_adj = (S - divPV) * Z_csa * (1 + 0.01 * (fwd_rate(dt_stock_settle, dt_exp_settle, yc) + r_repo_sprd_pct) * (dt_exp_settle - dt_stock_settle).value / 360)
    σ = σ * S / (S - divPV)
    r_sofr2exp = log(df(dt_trade, yc) / df(dt_expiry, yc)) / ((dt_expiry - dt_trade).value / 365)
    r_csa2exp = r_sofr2exp + 0.01 * r_csa_sprd_pct
    d1 = (log(S_adj / K) + (r_csa2exp + σ^2 / 2) * t) / (σ * √t)
    #  = (log(S*exp((r_repo-r_csa)*t)/K) + (r_csa + σ^2/2)*t)/(σ * √t)
    #   = (log(S/K)+(r_repo - r_csa)*t + (r_csa +σ^2/2)*t)/(σ * √t)
    #   = (log(S/K) + (r_repo + +σ^2/2)*t)/(σ * √t)
    d2 = d1 - σ * √t
    dist = Normal()
    if lowercase(optType)[1] == 'c'
        return S_adj * cdf(dist, d1) - K * cdf(dist, d2) * Z_csa
    elseif lowercase(optType)[1] == 'p'
        return K * cdf(dist, -d2) * Z_csa - S_adj * cdf(dist, -d1)
    end
end

function implied_vol_BS_yc(optPrice, yc, optType, S, K, r_csa_sprd_pct, r_repo_sprd_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
    obj(x) = (BS_yc(yc, optType, S, K, r_csa_sprd_pct, r_repo_sprd_pct, x, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt) - optPrice)^2
    res = Optim.optimize(obj, 0, 1000)
    vol = Optim.minimizer(res)
    vol
end

function conversion_BS_yc(yc, S, K, r_csa_sprd_pct, r_repo_sprd_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
    callPrice = BS_yc(yc, "c", S, K, r_csa_sprd_pct, r_repo_sprd_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
    putPrice = BS_yc(yc, "p", S, K, r_csa_sprd_pct, r_repo_sprd_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
    conversion = putPrice - callPrice + S - K
    return (conversion)
end

function implied_repo_spread_from_conversion_BS_yc(convPrice, yc, S, K, r_csa_sprd_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt)
    T = (dt_expiry - dt_trade).value

    #= the BS function creates the requisite continuos compounding repo rate as follows:
    days = (dt_expiry - dt_trade).value
    t = days / 365
    r_repo = log(1 + 0.01 * r_repo_mm_pct * days / 360) / t
    we want the argument of log to remain positive so we require
    1 + .01 * r_repo_mm_pct * days / 360 > 0
    r_repo_mm_pct > -360 * 100 / days
    =#

    (ycinfo, futs, intervals) = yc
    max_allowable_repo_sprd_pct = maximum(intervals.rate) + 5.0
    min_allowable_repo_sprd_pct = maximum(intervals.rate) - 36000 / T

    obj(x) = (conversion_BS_yc(yc, S, K, r_csa_sprd_pct, x, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt) - convPrice)^2
    res = Optim.optimize(obj, min_allowable_repo_sprd_pct, max_allowable_repo_sprd_pct)
    if Optim.converged(res)
        return Optim.minimizer(res)
    else
        return NaN
    end
end

function crr_yc(yc, optType, S, K, r_csa_sprd_pct, r_repo_sprd_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)

    (ycinfo, futs, intervals) = yc
    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end
    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)
    dt_exp_settle = workday(dt_expiry, exercise_settlement)

    # put in error trap code in case ycinfo.date[end] < dt_exp_settle

    ndays = (dt_expiry - dt_trade).value
    ndays2 = (dt_exp_settle - dt_trade).value
    T = ndays / 360
    N = ceil(Int, stepsPerDay * ndays) + 2
    N2 = ceil(Int, stepsPerDay * ndays2) + 2
    Δt = T / (N - 2)  # N.B. already in units of 360 day year

    rateCC = 360 * 100 * log.(1 .+ .01 * ycinfo.rate ./ 360)

    dateVec = dt_trade .+ Day.(ceil.(round.((-2:N2-2) * 360 * Δt, digits=1)))
    dateVecNext = vcat(dateVec[2:end], dateVec[end])
    i_stock_settle = maximum(findall(dateVec .== dt_stock_settle))
    i_option_settle = maximum(findall(dateVec .== dt_option_settle))
    baseRateVec = zeros(length(dateVec)) .+ NaN
    dates = unique(dateVec[dateVec.>=dt_trade])
    for d ∈ dates
        rng1 = findall(dateVec .== d)
        baseRateVec[rng1] .= rateCC[findall(ycinfo.date .== d)[1]]
    end
    stIndx = findall(dateVec .== dt_trade)[1]
    csaVec = baseRateVec .+ r_csa_sprd_pct
    # so that we only discount back to the option settlement date
    csaVec[findall(dateVec .<= dt_option_settle)] .= 0
    repoVec = baseRateVec .+ r_repo_sprd_pct
    # so that repo accrual begins on stock_settle date
    repoVec[findall(dateVec .<= dt_stock_settle)] .= 0
    dfCsaVec = ones(length(dateVec))
    dfCsaVec[stIndx:end] = [1; exp.(-cumsum(0.01 * csaVec[stIndx:end-1] * Δt))]
    dfRepoVec = ones(length(dateVec))
    dfRepoVec[stIndx:end] = [1; exp.(-cumsum(0.01 * repoVec[stIndx:end-1] * Δt))]
    aVec = exp.(0.01 * repoVec * Δt)
    Rvec = exp.(0.01 * csaVec * Δt)
    # so that we only discount back to the option settlement date
    Rvec[findall(dateVecNext .<= dt_option_settle)] .= 1

    σ = vol_pct / 100
    divfvs = zeros(Float64, N + 1)
    if sum(divamt) > 0
        # we want at each slice divfv = value of all divs w/ exdiv not in the past as of that slice
        i_exdiv = zeros(Int, length(divamt))
        i_divpmt = zeros(Int, length(divamt))
        for i = 1:length(divamt)
            i_divpmt[i] = floor(Int, 4 + stepsPerDay * ((dt_divpmt[i] - dt_trade).value - 1))
            i_exdiv[i] = floor(Int, 4 + stepsPerDay * ((dt_exdiv[i] - dt_trade).value - 1))
            # this gets discounting right even when divpmt falls between slices because stepsPerDay < 1
            #divfvs[1:i_exdiv[i]-1] .+= (divamt[i] * dfRepoVec[i_divpmt[i]]) ./ dfRepoVec[1:i_exdiv[i]-1]
            divfvs[1:i_exdiv[i]-1] .+= (divamt[i] * df(dt_divpmt[i], yc)) ./ dfCsaVec[1:i_exdiv[i]-1]
        end
    end
    # we will discount using R_csa from i_option_settle to end
    # we will grow S0 at R_repo from i_stock_settle to end
    # this means we'll go to the LAST slice on dt_option_settle and dt_stock_settle
    #i_option_settle = 1 + stepsPerDay * (dt_option_settle - dt_trade).value
    #i_stock_settle = 1 + stepsPerDay * (dt_stock_settle - dt_trade).value
    S0 = S - divfvs[1]
    σ = σ * S / S0

    U = exp(σ * √Δt)
    D = 1.0 / U
    pVec = (aVec .- D) / (U - D)
    # so that we only start repo accrual after stock settle
    pVec[findall(dateVecNext .<= dt_stock_settle)] .= (1 - D) / (U - D)
    qVec = 1 .- pVec

    local Z_csa::Float64
    local Z_repo::Float64
    local V::Vector{Float64}
    local p2use::Float64
    local q2use::Float64
    local R2use::Float64
    local px::Float64
    local x::Float64
    local y::Float64
    local tau::Float64

    if lowercase(optType)[1] == 'p'
        dt = dt_expiry
        dt_settle = workday(dt, exercise_settlement)
        i_settle = maximum(findall(dateVec .== dt_settle))
        Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, N + 1)]
        Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, N + 1)]
        #Z_csa = dfCsaVec[i_settle] / dfCsaVec[N+1]
        #Z_repo = dfRepoVec[i_settle] / dfRepoVec[N+1]

        V = Z_csa * [max(0.0, K - S0 * U^i * D^(N - i) / Z_repo) for i = 0:N]
        for n = N-1:-1:0
            dt = dt_trade + Day(ceil((n - 2) / (N - 2) * ndays))
            dt_settle = workday(dt, exercise_settlement)
            i_settle = maximum(findall(dateVec .== dt_settle))
            Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, n + 1)]
            Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, n + 1)]
            #Z_csa = dfCsaVec[i_settle] / dfCsaVec[n+1]
            #Z_repo = dfRepoVec[i_settle] / dfRepoVec[n+1]

            if dt == bd(dt)
                for i = 0:n
                    y = (qVec[n+1] * V[i+1] + pVec[n+1] * V[i+2]) / Rvec[n+1]
                    px = S0 * U^i * D^(n - i) + divfvs[n+1]
                    x = Z_csa * (K - px / Z_repo)
                    V[i+1] = max(x, y)
                end
            else
                for i = 0:n
                    V[i+1] = (qVec[n+1] * V[i+1] + pVec[n+1] * V[i+2]) / Rvec[n+1]
                end
            end
            if n == 2
                delta = (V[3] - V[1]) / (S0 * U^2 - S0 * D^2)
                gamma = 2 * ((V[3] - V[2]) / (S0 * U^2 - S0 * U * D) - (V[2] - V[1]) / (S0 * U * D - S0 * D^2)) / (S0 * U^2 - S0 * D^2)
                value = V[2]
                return (value, delta, gamma)
            end
        end
    elseif lowercase(optType)[1] == 'c'
        dt = dt_expiry
        dt_settle = workday(dt, exercise_settlement)
        i_settle = maximum(findall(dateVec .== dt_settle))
        Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, N + 1)]
        Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, N + 1)]
        #Z_csa = dfCsaVec[i_settle] / dfCsaVec[N+1]
        #Z_repo = dfRepoVec[i_settle] / dfRepoVec[N+1]
        V = Z_csa * [max(0.0, S0 * U^i * D^(N - i) / Z_repo - K) for i = 0:N]
        for n = N-1:-1:0
            dt = dt_trade + Day(ceil((n - 2) / (N - 2) * ndays))
            dt_settle = workday(dt, exercise_settlement)
            i_settle = maximum(findall(dateVec .== dt_settle))
            Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, n + 1)]
            Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, n + 1)]
            #Z_csa = dfCsaVec[i_settle] / dfCsaVec[n+1]
            #Z_repo = dfRepoVec[i_settle] / dfRepoVec[n+1]
            if dt == bd(dt)
                for i = 0:n
                    y = (qVec[n+1] * V[i+1] + pVec[n+1] * V[i+2]) / Rvec[n+1]
                    px = S0 * U^i * D^(n - i) + divfvs[n+1]
                    x = Z_csa * (px / Z_repo - K)
                    V[i+1] = max(x, y)
                end
            else
                for i = 0:n
                    V[i+1] = (qVec[n+1] * V[i+1] + pVec[n+1] * V[i+2]) / Rvec[n+1]
                end
            end
            if n == 2
                delta = (V[3] - V[1]) / (S0 * U^2 - S0 * D^2)
                gamma = 2 * ((V[3] - V[2]) / (S0 * U^2 - S0 * U * D) - (V[2] - V[1]) / (S0 * U * D - S0 * D^2)) / (S0 * U^2 - S0 * D^2)
                value = V[2]
                return (value, delta, gamma)
            end
        end
    end
end

function trinomial_yc(yc, optType, S, K, r_csa_sprd_pct, r_repo_sprd_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)

    (ycinfo, futs, intervals) = yc
    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end
    #exercise_settlement = 2
    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)
    dt_exp_settle = workday(dt_expiry, exercise_settlement)

    # put in error trap code in case ycinfo.date[end] < dt_exp_settle

    ndays = (dt_expiry - dt_trade).value
    ndays2 = (dt_exp_settle - dt_trade).value
    T = ndays / 360
    N = ceil(Int, stepsPerDay * ndays) + 1
    N2 = ceil(Int, stepsPerDay * ndays2) + 1
    Δt = T / (N - 1)  # N.B. already in units of 360 day year

    rateCC = 360 * 100 * log.(1 .+ 0.01 * ycinfo.rate ./ 360)

    dateVec = dt_trade .+ Day.(ceil.(round.((-1:N2-1) * 360 * Δt, digits=1)))
    dateVecNext = vcat(dateVec[2:end], dateVec[end])
    i_stock_settle = maximum(findall(dateVec .== dt_stock_settle))
    i_option_settle = maximum(findall(dateVec .== dt_option_settle))
    baseRateVec = zeros(length(dateVec)) .+ NaN
    dates = unique(dateVec[dateVec.>=dt_trade])
    for d ∈ dates
        rng1 = findall(dateVec .== d)
        baseRateVec[rng1] .= rateCC[findall(ycinfo.date .== d)[1]]
    end
    stIndx = findall(dateVec .== dt_trade)[1]
    csaVec = baseRateVec .+ r_csa_sprd_pct
    # so that we only discount back to the option settlement date
    csaVec[findall(dateVecNext .<= dt_option_settle)] .= 0
    repoVec = baseRateVec .+ r_repo_sprd_pct
    # so that repo accrual begins on stock_settle date
    repoVec[findall(dateVecNext .<= dt_stock_settle)] .= 0
    dfCsaVec = ones(length(dateVec))
    dfCsaVec[stIndx:end] = [1; exp.(-cumsum(0.01 * csaVec[stIndx:end-1] * Δt))]
    dfRepoVec = ones(length(dateVec))
    dfRepoVec[stIndx:end] = [1; exp.(-cumsum(0.01 * repoVec[stIndx:end-1] * Δt))]
    Rvec = exp.(0.01 * csaVec * Δt)

    σ = vol_pct / 100
    divfvs = zeros(Float64, N + 1)
    if sum(divamt) > 0
        # we want at each slice divfv = value of all divs w/ exdiv not in the past as of that slice
        i_exdiv = zeros(Int, length(divamt))
        i_divpmt = zeros(Int, length(divamt))
        for i = 1:length(divamt)
            i_divpmt[i] = floor(Int, 3 + stepsPerDay * ((dt_divpmt[i] - dt_trade).value - 1))
            i_exdiv[i] = floor(Int, 3 + stepsPerDay * ((dt_exdiv[i] - dt_trade).value - 1))
            # this gets discounting right even when divpmt falls between slices because stepsPerDay < 1
            #divfvs[1:i_exdiv[i]-1] .+= (divamt[i] * dfRepoVec[i_divpmt[i]]) ./ dfRepoVec[1:i_exdiv[i]-1]
            divfvs[1:i_exdiv[i]-1] .+= (divamt[i] * df(dt_divpmt[i], yc)) ./ dfCsaVec[1:i_exdiv[i]-1]
        end
    end
    S0 = S - divfvs[1]
    σ = σ * S / S0

    u = exp(σ * sqrt(3 * Δt))
    m = 1.0
    d = 1.0 / u
    p_u_vec = 1 / 6 .+ (0.01 * repoVec .- σ^2 / 2) * sqrt(Δt / (12 * σ^2))

    local V::Vector{Float64}
    local px::Float64
    local x::Float64
    local y::Float64
    local Z_csa::Float64
    local Z_repo::Float64
    local i_settle::Int32

    if lowercase(optType)[1] == 'p'
        dt = dt_expiry
        dt_settle = workday(dt, exercise_settlement)
        i_settle = maximum(findall(dateVec .== dt_settle))
        i_settle = findall(dateVec .== dt_settle)[end]
        Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, N + 1)]
        Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, N + 1)]
        V = Z_csa * [max(0.0, K - S0 * u^i / Z_repo) for i = -N:N]
        for n = N-1:-1:0
            dt = dt_trade + Day(ceil((n - 1) / (N - 1) * ndays))
            dt_settle = workday(dt, exercise_settlement)
            i_settle = maximum(findall(dateVec .== dt_settle))
            i_settle = findall(dateVec .== dt_settle)[end]
            Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, n + 1)]
            Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, n + 1)]
            if dt == bd(dt)
                for i = -n:n
                    y = ((1 - p_u_vec[n+1] - 2 / 3) * V[i+n+1] + 2 / 3 * V[i+n+2] + p_u_vec[n+1] * V[i+n+3]) / Rvec[n+1]
                    px = S0 * u^i + divfvs[n+1]
                    x = Z_csa * (K - px / Z_repo)
                    V[i+n+1] = max(x, y)
                end
            else
                for i = -n:n
                    V[i+n+1] = ((1 - p_u_vec[n+1] - 2 / 3) * V[i+n+1] + 2 / 3 * V[i+n+2] + p_u_vec[n+1] * V[i+n+3]) / Rvec[n+1]
                end
            end
            if n == 1
                delta = (V[3] - V[1]) / (S0 * u - S0 * d)
                gamma = 2 * ((V[3] - V[2]) / (S0 * u - S0 * m) - (V[2] - V[1]) / (S0 * m - S0 * d)) / (S0 * u - S0 * d)
                value = V[2]
                return (value, delta, gamma)
            end
        end
    elseif lowercase(optType)[1] == 'c'
        dt = dt_expiry
        dt_settle = workday(dt, exercise_settlement)
        i_settle = maximum(findall(dateVec .== dt_settle))
        i_settle = findall(dateVec .== dt_settle)[end]
        Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, N + 1)]
        Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, N + 1)]
        V = Z_csa * [max(0.0, S0 * u^i / Z_repo - K) for i = -N:N]
        for n = N-1:-1:0
            dt = dt_trade + Day(ceil((n - 1) / (N - 1) * ndays))
            dt_settle = workday(dt, exercise_settlement)
            i_settle = maximum(findall(dateVec .== dt_settle))
            i_settle = findall(dateVec .== dt_settle)[end]
            Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, n + 1)]
            Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, n + 1)]
            if dt == bd(dt)
                for i = -n:n
                    y = ((1 - p_u_vec[n+1] - 2 / 3) * V[i+n+1] + 2 / 3 * V[i+n+2] + p_u_vec[n+1] * V[i+n+3]) / Rvec[n+1]
                    px = S0 * u^i + divfvs[n+1]
                    x = Z_csa * (px / Z_repo - K)
                    V[i+n+1] = max(x, y)
                end
            else
                for i = -n:n
                    V[i+n+1] = ((1 - p_u_vec[n+1] - 2 / 3) * V[i+n+1] + 2 / 3 * V[i+n+2] + p_u_vec[n+1] * V[i+n+3]) / Rvec[n+1]
                end
            end
            if n == 1
                delta = (V[3] - V[1]) / (S0 * u - S0 * d)
                gamma = 2 * ((V[3] - V[2]) / (S0 * u - S0 * m) - (V[2] - V[1]) / (S0 * m - S0 * d)) / (S0 * u - S0 * d)
                value = V[2]
                return (value, delta, gamma)
            end
        end
    end
end

function option_am_yc(yc, treeType, optType, S, K, r_csa_sprd_pct, r_repo_sprd_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    if lowercase(treeType) == "crr"
        fn = crr_yc
    elseif lowercase(treeType) == "trinomial"
        fn = trinomial_yc
    end
    return fn(yc, optType, S, K, r_csa_sprd_pct, r_repo_sprd_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)#[1]
end

function implied_vol_am_yc(optPrice, yc, treeType, optType, S, K, r_csa_sprd_pct, r_repo_sprd_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    obj(x) = (option_am_yc(yc, treeType, optType, S, K, r_csa_sprd_pct, r_repo_sprd_pct, x, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1] - optPrice)^2
    res = Optim.optimize(obj, 0, 1000)
    vol = Optim.minimizer(res)
    vol
end

function implied_vol_and_repo_from_put_and_call_am_yc(yc, treeType, putPrice, callPrice, S, K, r_csa_sprd_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    # this function simultaneously solves for the implied volatility and repo spread for a pair of same-strike options at a given csa spread
    # VERY USEFUL if you don't have a good guess for the repo rate from your Prime Broker

    obj(x) = (callPrice - option_am_yc(yc, treeType, "call", S, K, r_csa_sprd_pct, x[1], x[2]^2, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1])^2 +
             (putPrice - option_am_yc(yc, treeType, "put", S, K, r_csa_sprd_pct, x[1], x[2]^2, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1])^2

    initial_x = [r_csa_sprd_pct, sqrt(50)]
    res = Optim.optimize(obj, initial_x)
    if Optim.converged(res)
        r_repo_mm_pct = Optim.minimizer(res)[1]
        vol_pct = Optim.minimizer(res)[2]^2
        return (r_repo_mm_pct, vol_pct)
        #println("call and put are both correctly priced when:")
        #println("r_repo_mm_pct = $(round(r_repo_mm_pct,digits=2)) and")
        #println("vol_pct = $(round(vol_pct))\n")
    else
        return (NaN, NaN)
        #println("Optimizer did not find a solution")
    end
end

function conversion_am_yc(yc, treeType, S, K, r_csa_sprd_pct, r_repo_sprd_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    call = option_am_yc(yc, treeType, "c", S, K, r_csa_sprd_pct, r_repo_sprd_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    put = option_am_yc(yc, treeType, "p", S, K, r_csa_sprd_pct, r_repo_sprd_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    value = put[1] - call[1] + S - K
    delta = put[2] - call[2] + 1
    gamma = put[3] - call[3]
    return (value, delta, gamma)
end


function implied_dividend_from_conversion_am_yc(convPrice, yc, treeType, S, K, r_csa_sprd_pct, r_repo_sprd_pct, vol_pct, dt_trade, dt_expiry, dt_exdiv, dt_divpmt, stepsPerDay, maxDiv2Px)
    n = length(dt_exdiv)
    obj(x) = (conversion_am_yc(yc, treeType, S, K, r_csa_sprd_pct, r_repo_sprd_pct, vol_pct, dt_trade, dt_expiry, x * ones(n), dt_exdiv[1:n], dt_divpmt[1:n], stepsPerDay)[1] - convPrice)^2
    ans = NaN
    maxrng = 0.01
    while (maxrng < maxDiv2Px)
        res = Optim.optimize(obj, S * (maxrng - 0.01), S * maxrng, abs_tol=0.001)
        if Optim.converged(res) && Optim.minimum(res) < 0.0001
            ans = Optim.minimizer(res)
            break
        else
            maxrng = maxrng + 0.01
        end
    end
    return ans
end

function implied_repo_spread_from_conversion_am_yc(convPrice, yc, treeType, S, K, r_csa_sprd_pct, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)

    #=
    the put/call valuations use r_repo = log(1.0 + 0.01 * r_repo_mm_pct * ndays / 360) / T
    so we need to require r_repo_mm_pct > 1 + -36000/ndays
    =#

    T = (dt_expiry - dt_trade).value
    (ycinfo, futs, intervals) = yc
    global max_repo_sprd_pct = maximum(intervals.rate) + 5.0
    global min_repo_sprd_pct = maximum(intervals.rate) - 36000 / T

    offset_vec = vcat(10, 0, -1:-1:-3, -10:-10:-100, -200:-100:-900, -1000:-1000:-3000)
    offset_vec = min.(offset_vec, max_repo_sprd_pct)
    offset_vec = max.(offset_vec, min_repo_sprd_pct)
    offset_vec = sort(unique(offset_vec), rev=true)
    pvec = zeros(length(offset_vec))
    for i ∈ eachindex(offset_vec)
        pvec[i] = conversion_am_yc(yc, treeType, S, K, r_csa_sprd_pct, offset_vec[i], vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1]
        if (i > 1) && ((pvec[i] - convPrice) * (pvec[i-1] - convPrice) <= 0) # i.e. the soln lies in between rvec[i] and rvec[i-1]
            max_repo_sprd_pct = offset_vec[i-1]
            min_repo_sprd_pct = offset_vec[i]
            break
        elseif (pvec[i] > convPrice)
            return NaN # pvec is a decreasing fn of r so as we walk r in the negative direction, pvec > mktPrice means we can't cross from below and should exit
        end
    end
    if abs(pvec[end]) > 0
        return NaN # pvec[end] will be non-zero only if no soln was found
    end
    obj(x) = (conversion_am_yc(yc, treeType, S, K, r_csa_sprd_pct, x, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1] - convPrice)^2
    res = Optim.optimize(obj, min_repo_sprd_pct, max_repo_sprd_pct)
    if Optim.converged(res)
        return Optim.minimizer(res)
    end
end

function implied_repo_spread_from_jelly_roll_yc(mktPrice, yc, treeType, S, r_csa_sprd_pct, dt_trade, divamt, dt_exdiv, dt_divpmt, stepsPerDay, K1, vol_pct1, dt_expiry1, K2, vol_pct2, dt_expiry2, max_repo_sprd_pct=10)
    (ycinfo, futs, intervals) = yc
    T1 = (dt_expiry1 - dt_trade).value
    T2 = (dt_expiry2 - dt_trade).value
    min_repo_sprd_pct = maximum(yc[3].rate) - 36000 / max(T1, T2)

    offset_vec =
        vcat(
            max_repo_sprd_pct,
            -1:-1:-3,
            -10:-10:-100,
            -200:-100:-900, -1000:-1000:-3000
        )
    offset_vec = offset_vec[offset_vec.>=min_repo_sprd_pct]
    pvec = zeros(length(offset_vec))
    for i ∈ eachindex(offset_vec)
        pvec[i] =
            conversion_am_yc(yc, treeType, S, K2, r_csa_sprd_pct, offset_vec[i], vol_pct2, dt_trade, dt_expiry2, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1] -
            conversion_am_yc(yc, treeType, S, K1, r_csa_sprd_pct, offset_vec[i], vol_pct1, dt_trade, dt_expiry1, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1]
    end
    global max_repo_sprd = []
    global min_repo_sprd = []
    for i = 2:length(offset_vec)
        if (pvec[i] - mktPrice) * (pvec[i-1] - mktPrice) <= 0 # i.e. the soln lies in between rvec[i] and rvec[i-1]
            global max_repo_sprd = [max_repo_sprd; offset_vec[i-1]]
            global min_repo_sprd = [min_repo_sprd; offset_vec[i]]
        end
    end
    global possAns = []
    obj(x) = (conversion_am_yc(yc, treeType, S, K2, r_csa_sprd_pct, x, vol_pct2, dt_trade, dt_expiry2, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1] -
              conversion_am_yc(yc, treeType, S, K1, r_csa_sprd_pct, x, vol_pct1, dt_trade, dt_expiry1, divamt, dt_exdiv, dt_divpmt, stepsPerDay)[1] -
              mktPrice)^2
    if length(max_repo_sprd) > 0
        for i ∈ eachindex(max_repo_sprd)
            res = Optim.optimize(
                obj,
                min_repo_sprd[i],
                max_repo_sprd[i],
                abs_tol=1e-5,
            )
            if Optim.converged(res)
                global possAns = [possAns; Optim.minimizer(res)]
            end
        end
    end
    if length(possAns) == 0
        return NaN
    else
        return possAns #[findmin(identity, abs.(possAns .- r_csa))[2]]
    end
end

function yc_hedger(f, yc, baseRate="SOFR", dt_trade=today(), dEnd=today() + Year(2), bumpSize=0.1)
    #info = curveStripperFromFuts(baseRate, dt_trade, dEnd)
    (ycinfo, futs, intervals) = yc
    nf = nrow(futs)
    val = f(yc)

    vals = vcat(val, zeros(nf) .+ NaN)
    hvals = zeros(nf, nf + 1) .+ NaN

    # first col of hvals will be fut prices evaluated with base curve
    # following cols will be fut prices evaluated with perturbed curves

    # vals[1] will contain caplet val from base curve
    # vals [i+1] will contain the caplet val obtained by bumping by 1bp
    # the instantaneous fwds in info between intervals.startDate[i] and intervals.endDate[i]
    for i = 1:nf+1
        global yc = copy(ycinfo)
        if i > 1
            if i == 2
                d_start = dt_trade
            else
                d_start = futs.lastRateDate[i-2] + Day(1)
            end
            d_finish = futs.lastRateDate[i-1]
            rng = (findall(yc.date .>= d_start .&& yc.date .<= d_finish))
            yc.rate[rng] = yc.rate[rng] .+ bumpSize
            # recompute discount factors in yc now that we've bumped rates
            r = 2
            while r <= nrow(yc)
                if !isnan(yc.df[r])
                    r_start = max(1, r - 5)
                    r_end = r - 1
                    r_prev = r_start - 1 + maximum(findall(yc.df[r_start:r_end] .> 0))
                    yc.df[r] = yc.df[r_prev] / (1 + 0.01 * yc.rate[r_prev] * (yc.date[r] - yc.date[r_prev]).value / 360)
                end
                r = r + 1
            end
        end
        vals[i] = f((yc, futs, intervals))
        for r = 1:nrow(futs)
            hvals[r, i] = impliedFutPx(string(futs.bct[r]), (yc, futs, intervals))
        end
    end
    delVals = (vals[2:end] .- vals[1])
    delHvals = hvals[:, 2:end] .- hvals[:, 1]
    pointVals = zeros(nf)
    for r = 1:nf
        if futs.bct[r][1:2] == "sq"
            pointVals[r] = 2500
        elseif futs.bct[r][1:2] ∈ ["sl", "zq"]
            pointVals[r] = 4167
        end
        delHvals[r, :] = delHvals[r, :] .* pointVals[r]
    end
    delHvals = delHvals'
    delVals[abs.(delVals).<1e-3] .= 0
    delHvals[abs.(delHvals).<1e-3] .= 0
    delVals = delVals ./ (bumpSize / 0.01)
    delHvals = delHvals ./ (bumpSize / 0.01)
    F = svd(delHvals)
    S = F.S
    S[S.==0] .= Inf
    hedges = F.Vt' * Diagonal(1 ./ S) * F.U' * (-delVals)
    return (round(val), round(sum(hedges .* pointVals) / 100), DataFrame(contract=futs.bct, hedge=round.(hedges, digits=1)))
    if false
        futs.hedgepv01 = hedges .* pointVals
        dates = d:Day(1):futs.lastRateDate[end]
        v = zeros(length(dates))
        for i = 1:nf
            indx = findall(dates .>= futs.firstRateDate[i] .&& dates .<= futs.lastRateDate[i])
            v[indx] = v[indx] .+ pointVals[i] * hedges[i]
        end
        plot(dates, v, ticks=:native)
    end
end

################ simple functions for benchmarking and illustration
function BSimple(optType::String, S::Union{Float64,Int}, K::Union{Float64,Int}, μ::Float64, r::Float64, q::Union{Float64,Int}, σ_bs::Float64, T::Union{Float64,Int}, rebalsPerDay::Int=1)
    Δt = 1 / 365 / rebalsPerDay
    #σ_adjusted = σ * (1 + Δt * (μ - r) * (r - μ - σ^2) / (2 * σ^2))
    d1 = (log(S / K) + (r - q + σ_bs^2 / 2) * T) / (σ_bs * √T)
    d2 = d1 - σ_bs * √T
    dist = Normal()
    if lowercase(optType)[1] == 'c'
        V = exp(-q * T) * S * cdf(dist, d1) - K * exp(-r * T) * cdf(dist, d2)
        Δ = exp(-q * T) * cdf(dist, d1)
        Γ = exp(-q * T) * pdf(dist, d1) / (S * σ_bs * √T)
        Θ = -S * pdf(dist, d1) * σ_bs / (2√T) - r * K * exp(-(r - q) * T) * cdf(dist, d2)
        vega = S * √T * pdf(dist, d1)
    elseif lowercase(optType)[1] == 'p'
        V = exp(-r * T) * K * cdf(dist, -d2) - S * exp(-q * T) * cdf(dist, -d1)
        Δ = exp(-q * T) * (cdf(dist, d1) - 1)
        Γ = exp(-q * T) * pdf(dist, d1) / (S * σ_bs * √T)
        Θ = -S * pdf(dist, d1) * σ_bs / (2√T) + r * K * exp(-(r - q) * T) * cdf(dist, -d2)
        vega = S * √T * pdf(dist, d1)
    end
    Δ_adjusted = Δ + (μ - r + σ_bs^2 / 2) * Δt * S * Γ
    return (V, Δ, Γ, Θ, vega, Δ_adjusted)
end

function B76(optType, F, K, r, σ_bs, T)
    d1 = (log(F / K) + ( σ_bs^2 / 2) * T) / (σ_bs * √T)
    d2 = d1 - σ_bs * √T
    dist = Normal()
    if lowercase(optType)[1] == 'c'
        V = exp(-r * T) * (F * cdf(dist, d1) - K * cdf(dist, d2))
        Δ = exp(-r * T) * cdf(dist, d1)
        Γ = exp(-r * T) * pdf(dist, d1) / (F * σ_bs * √T)
        Θ = -F * pdf(dist, d1) * σ_bs / (2√T) - r * K * cdf(dist, d2)
        vega = F * √T * pdf(dist, d1)
    elseif lowercase(optType)[1] == 'p'
        V = exp(-r * T) * (K * cdf(dist, -d2) - F * cdf(dist, -d1))
        Δ = exp(-r * T) * (cdf(dist, d1) - 1)
        Γ = exp(-r * T) * pdf(dist, d1) / (F * σ_bs * √T)
        Θ = -F * pdf(dist, d1) * σ_bs / (2√T) + r * K * cdf(dist, -d2)
        vega = F * √T * pdf(dist, d1)
    end
    return (V, Δ, Γ, Θ, vega)
end

function implied_vol_B76(optPrice, optType, F, K, r, T)
    obj(x) = (B76(optType, F, K,r, x, T)[1] - optPrice)^2
    res = Optim.optimize(obj, 0, 20)
    vol = Optim.minimizer(res)
    vol
end

function implied_vol_BSimple(optPrice, optType, S, K, μ, r, q, T, rebalsPerDay)
    obj(x) = (BSimple(optType, S, K, μ, r, q, x, T, rebalsPerDay)[1] - optPrice)^2
    res = Optim.optimize(obj, 0, 20)
    vol = Optim.minimizer(res)
    vol
end

#=
function conversion_BS_simple(S, K, μ, r, q, σ, T, rebalsPerDay)
    callPrice = BSimple("c", S, K, μ, r, q, σ, T, rebalsPerDay)[1]
    putPrice = BSimple("p", S, K, μ, r, q, σ, T, rebalsPerDay)[1]
    conversion = putPrice - callPrice + S - K
    return (conversion)
end

function implied_repo_from_conversion_BS_simple(convPrice, S, K, μ, r, σ, T, rebalsPerDay)
    global max_allowable_q = r  #recall q = r - r_repo
    global min_allowable_q = r

    offset_vec = vcat(10, 0, -1:-1:-3, -10:-10:-100, -200:-100:-900, -1000:-1000:-3000) ./ 100
    pvec = zeros(length(offset_vec))
    qvec = r .- offset_vec
    for i = 1:length(qvec)
        pvec[i] = conversion_BS_simple(S, K, μ, r, qvec[i], σ, T, rebalsPerDay)
        if (i > 1) && ((pvec[i] - convPrice) * (pvec[i-1] - convPrice) < 0) # i.e. the soln lies in between rvec[i] and rvec[i-1]
            min_allowable_q = qvec[i-1]
            max_allowable_q = qvec[i]
            break
        elseif (pvec[i] > convPrice)
            return NaN # pvec is a decreasing fn of r so as we walk r in the negative direction, pvec > mktPrice means we can't cross from below and should exit
        end
    end
    if abs(pvec[end]) > 0
        return NaN # pvec[end] will be non-zero only if no soln was found
    end
    obj(x) = (conversion_BS_simple(S, K, μ, r, x, σ, T, rebalsPerDay) - convPrice)^2
    res = Optim.optimize(obj, min_allowable_q, max_allowable_q)
    if Optim.converged(res)
        return 100 * (r - Optim.minimizer(res))
    end
end
=#

function cev_simple(optType::String, S::Union{Float64,Int}, K::Union{Float64,Int}, μ::Float64, r::Float64, q::Union{Float64,Int}, σ_bs, T::Union{Float64,Int}, β::Float64, rebalsPerDay::Int=1)
    Δt = 1 / 365 / rebalsPerDay
    #if β > 1
    #    return ("ERROR:  function only works for β < 1")
    #end
    # adopting the notation of the Larguinho et al paper but using K rather that X and β rather than β/2
    F(w, ν, λ) = cdf(NoncentralChisq(ν, λ), w)
    Q(w, ν, λ) = 1 - F(w, ν, λ)
    p(w, ν, λ) = pdf(NoncentralChisq(ν, λ), w)
    σ = σ_bs * S^(1 - β)
    k = (r - q) / (σ^2 * (1 - β)) / (exp((r - q) * (2 - 2β) * T) - 1)
    x = k * S^(2 - 2β) * exp((r - q) * (2 - 2β) * T)
    y = k * K^(2 - 2β)

    if lowercase(optType)[1] == 'c'
        if β < 1
            V = S * exp(-q * T) * Q(2y, 2 + 1 / (1 - β), 2x) - K * exp(-r * T) * F(2x, 1 / (1 - β), 2y)
            Δ = exp(-q * T) * Q(2y, 2 + 1 / (1 - β), 2x) + 4x * (1 - β) / S * (S * exp(-q * T) * p(2y, 4 + 1 / (1 - β), 2x) - K * exp(-r * T) * p(2x, 1 / (1 - β), 2y))
            Θ = S * q * exp(-q * T) * Q(2y, 2 + 1 / (1 - β), 2x) - K * r * exp(-r * T) * F(2x, 1 / (1 - β), 2y) + 4x * (r - q) * (1 - β) / (-1 + exp((r - q) * (2 - 2β) * T)) * (S * exp(-q * T) * p(2y, 4 + 1 / (1 - β), 2x) - K * exp(-r * T) * p(2x, 1 / (1 - β), 2y))
        elseif β > 1
            V = S * exp(-q * T) * Q(2x, 1 / (β - 1), 2y) - K * exp(-r * T) * F(2y, 2 + 1 / (β - 1), 2x)
            Δ = exp(-q * T) * Q(2x, 1 / (β - 1), 2y) + 4x * (β - 1) / S * (S * exp(-q * T) * p(2x, 1 / (β - 1), 2y) - K * exp(-r * T) * p(2y, 4 + 1 / (β - 1), 2x))
            Θ = S * q * exp(-q * T) * Q(2x, 1 / (β - 1), 2y) - K * r * exp(-r * T) * F(2y, 2 + 1 / (β - 1), 2x) + 4x * (r - q) * (β - 1) / (-1 + exp((r - q) * (2 - 2β) * T)) * (S * exp(-q * T) * p(2x, 1 / (β - 1), 2y) - K * exp(-r * T) * p(2y, 4 + 1 / (β - 1), 2x))
        end
    elseif lowercase(optType)[1] == 'p'
        if β < 1
            V = exp(-r * T) * K - exp(-q * T) * S + S * exp(-q * T) * Q(2y, 2 + 1 / (1 - β), 2x) - K * exp(-r * T) * F(2x, 1 / (1 - β), 2y)
            Δ = -exp(-q * T) + exp(-q * T) * Q(2y, 2 + 1 / (1 - β), 2x) + 4x * (1 - β) / S * (S * exp(-q * T) * p(2y, 4 + 1 / (1 - β), 2x) - K * exp(-r * T) * p(2x, 1 / (1 - β), 2y))
            Θ = K * r * exp(-r * T) * Q(2x, 1 / (1 - β), 2y) - S * q * exp(-q * T) * F(2y, 2 + 1 / (1 - β), 2x) + 4x * (r - q) * (1 - β) / (-1 + exp((r - q) * (2 - 2β) * T)) * (S * exp(-q * T) * p(2y, 4 + 1 / (1 - β), 2x) - K * exp(-r * T) * p(2x, 1 / (1 - β), 2y))
        elseif β > 1
            V = exp(-r * T) * K - exp(-q * T) * S + S * exp(-q * T) * Q(2x, 1 / (β - 1), 2y) - K * exp(-r * T) * F(2y, 2 + 1 / (β - 1), 2x)
            Δ = -exp(-q * T) * F(2x, 1 / (β - 1), 2y) + 4x * (β - 1) / S * (S * exp(-q * T) * p(2x, 1 / (β - 1), 2y) - K * exp(-r * T) * p(2y, 4 + 1 / (β - 1), 2x))
            Θ = K * r * exp(-r * T) * Q(2y, 2 + 1 / (β - 1), 2x) - S * q * exp(-q * T) * F(2x, 1 / (β - 1), 2y) + 4x * (r - q) * (β - 1) / (-1 + exp((r - q) * (2 - 2β) * T)) * (S * exp(-q * T) * p(2x, 1 / (β - 1), 2y) - K * exp(-r * T) * p(2y, 4 + 1 / (β - 1), 2x))
        end
    end
    if β < 1
        Γ = 8x * (1 - β)^2 / S * exp(-q * T) * (((3 - 2β) / (2 - 2β) - x) * p(2y, 4 + 1 / (1 - β), 2x) + x * p(2y, 6 + 1 / (1 - β), 2x)) + 8x * (1 - β)^2 / S * exp(-r * T) * K / S * (x * p(2x, 1 / (1 - β), 2y) - y * p(2x, 2 + 1 / (1 - β), 2y))
        vega_bs = -4x / σ_bs * (S * exp(-q * T) * p(2y, 4 + 1 / (1 - β), 2x) - K * exp(-r * T) * p(2x, 1 / (1 - β), 2y))
        Δ_adjusted = Δ + (μ - r + σ_bs^2 / 2) * Δt * S
    elseif β > 1
        Γ = NaN
        vega_bs = NaN
        Δ_adjusted = Δ + (μ - r + σ_bs^2 / 2) * Δt * S
    end
    return (V, Δ, Γ, Θ, vega_bs, Δ_adjusted)
end

function crr_efficient_simple(optType::String, S::Union{Float64,Int}, K::Union{Float64,Int}, μ::Float64, r::Float64, q::Union{Float64,Int}, σ::Float64, T::Union{Float64,Int}, N::Int)
    #=
    let's first implement and check the extended binomial tree Pelsser and Vorst greek calculations described in Cruz-Dias paper
    as per figure 3 in that paper, we need to extend our binomial tree:  not only do we need the S(0,0) root and it's descendants
    but we need to start the tree at S(-2,-1) = S/(U*D) as the root and then build the tree from there so we can approximate delta by
    Δ = (V(0,1)) - V(0,-1)) / (S(0,1) -  S(0,-1)) and
    Γ = 2[(V(0,1) - V(0,0) / (S(0,1) - S(0,0)) - (V(0,0) - V(0,-1)) / (S(0,0) - S(0,-1))] / (S(0,1) - S(0,-1))
    =#
    if lowercase(optType)[1] == 'c'
        mult = 1
    else
        mult = -1
    end
    # we want sig * √T = sig'
    Δt = T / (N - 2)  # we need to extend the tree by 2 steps to get the right delta, but we don't want to increase the stdev on the terminal slice
    #Δt = T / N
    R = exp((r - q) * Δt)
    #σ = σ * (1 + Δt * (μ - r) * (r - μ - σ^2) / (2 * σ^2))
    U = exp(σ * √Δt)
    D = 1 / U
    p_up = (R - D) / (U - D)
    p_down = (U - R) / (U - D)

    V = zeros(N + 1)
    for i = 0:N
        V[i+1] = max(0, mult * (S * U^i * D^(N - i) - K))
    end
    thetaval1 = thetaval2 = Θ = value = Δ = Γ = Δ_adjusted = 0.0

    for n = N-1:-1:0
        for i = 0:n
            x = mult * (S * U^i * D^(n - i) - K)
            y = (p_down * V[i+1] + p_up * V[i+2]) / R
            V[i+1] = max(x, y)
        end
        if n == 4
            thetaval1 = V[3]
        end
        if n == 2
            value = V[2]
            Δ = (V[3] - V[1]) / (S * U^2 - S * D^2)
            Γ = 2 * ((V[3] - V[2]) / (S * U^2 - S * U * D) - (V[2] - V[1]) / (S * U * D - S * D^2)) / (S * U^2 - S * D^2)
        end
        if n == 0
            thetaval2 = V[1]
            Θ = (thetaval1 - thetaval2) / (4 * Δt)
            Δ_adjusted = Δ + (μ - r + σ^2 / 2) * Δt
            return (value, Δ, Γ, Θ, Δ_adjusted)
        end
    end
end

function trinomial_simple(optType::String, S::Union{Float64,Int}, K::Union{Float64,Int}, μ::Float64, r::Float64, q::Union{Float64,Int}, σ::Float64, T::Union{Float64,Int}, N::Int, λ=sqrt(2pi) / 2)
    #=
    we'll modify the Plesser-Vorst binomial tree logic for a trinomial tree
    =#
    if lowercase(optType)[1] == 'c'
        mult = 1
    else
        mult = -1
    end
    Δt = T / (N - 1)  # only need to extend trinomial tree by 1 step...
    R = exp(r * Δt)
    u = exp(λ * σ * sqrt(Δt))
    m = 1
    d = 1 / u
    r_repo = r - q
    M = exp(r_repo * Δt)
    Σ = M^2 * (-1 + exp(σ^2 * Δt))
    pᵤ = (u * (Σ + M^2 - M) - (M - 1)) / ((u - 1) * (u^2 - 1))
    p_d = (u^2 * (Σ + M^2 - M) - u^3 * (M - 1)) / ((u - 1) * (u^2 - 1))
    pₘ = 1 - pᵤ - p_d
    #u = exp(σ * sqrt(3 * Δt))
    #pᵤ = 1 / 6 + (r - q - σ^2 / 2) * sqrt(Δt / (12 * σ^2))
    #pₘ = 2 / 3
    #p_d  = 1-pₘ - pᵤ

    thetaval1 = thetaval2 = Θ = value = Δ = Γ = Δ_adjusted = 0.0
    V = zeros(2N + 1)
    for i = -N:N
        V[i+N+1] = max(0, mult * (S * u^i - K))
    end
    for n = N-1:-1:0
        for i = -n:n
            x = mult * (S * u^i - K)
            y = (p_d * V[i+n+1] + pₘ * V[i+n+2] + pᵤ * V[i+n+3]) / R
            V[i+n+1] = max(x, y)
        end
        if n == 2
            thetaval1 = V[3]
        end
        if n == 1
            Δ = (V[3] - V[1]) / (S * u - S * d)
            Γ = 2 * ((V[3] - V[2]) / (S * u - S * m) - (V[2] - V[1]) / (S * m - S * d)) / (S * u - S * d)
            value = V[2]
        end
        if n == 0
            thetaval2 = V[1]
            Θ = (thetaval1 - thetaval2) / (2 * Δt)
            Δ_adjusted = Δ + (μ - r + σ^2 / 2) * Δt
            return (value, Δ, Γ, Θ, Δ_adjusted)
        end
    end
end

function cev_efficient_simple(optType::String, s::Union{Float64,Int}, K::Union{Float64,Int}, μ::Float64, rt::Float64, q::Union{Float64,Int}, σ_bs::Float64, T::Union{Float64,Int}, N::Int, β::Float64)
    if lowercase(optType)[1] == 'c'
        mult = 1
    else
        mult = -1
    end
    Δt = T / (N - 2)  # we need to extend the tree by 2 steps to get the right delta, but we don't want to increase the stdev on the terminal slice
    R = exp(rt * Δt)
    #σ * S^β = σ_bs * S, so σ = σ_bs * S^(1-β)
    σ = σ_bs * s^(1 - β)
    α = 1 - β

    function x(S, σ, α)
        S^α / (α * σ)
    end

    function xInv(x, σ, α)
        # x = S^α / (α * σ)
        # so S = (x * α * σ)^(1/α)
        if x > 0
            return (x * α * σ)^(1 / α)
        else
            return 0
        end
    end

    function p_up(S, σ, α, rt, q, Δt)
        xval = x(S, σ, α)
        S_up = xInv(xval + √Δt, σ, α)
        S_down = xInv(xval - √Δt, σ, α)
        arg = (S * exp((rt - q) * Δt) - S_down) / (S_up - S_down)

        if xval > 0 && 0 <= arg <= 1
            return arg
        elseif xval <= 0 || arg < 0
            return 0
        elseif xval > 0 && arg > 1
            return 1
        end


    end

    s0 = s  ## this needs to be amended so that  S[3,2]=s in order to use the PV stuff
    x0 = x(s0, σ, α)
    # populate terminal values
    X = zeros(N + 1)
    S = 0 * X
    V = 0 * X
    thetaval1 = thetaval2 = Θ = value = Δ = Γ = Δ_adjusted = 0.0

    for i = 1:N+1
        X[i] = x0 + √Δt * (-N + 2 * (i - 1))
        S[i] = xInv(X[i], σ, α)
        V[i] = max(0, mult * (S[i] - K))
    end
    # propagate backward
    for n = N-1:-1:0
        for i = 0:n
            X[i+1] = x0 + √Δt * (-n + 2 * i)
            S[i+1] = xInv(X[i+1], σ, α)
            p = p_up(S[i+1], σ, α, rt, q, Δt)
            contVal = ((1 - p) * V[i+1] + p * V[i+2]) / R
            exVal = mult * (S[i+1] - K)
            V[i+1] = max(exVal, contVal)
        end
        if n == 4
            thetaval1 = V[3]
        end
        if n == 2
            value = V[2]
            Δ = (V[3] - V[1]) / (S[3] - S[1])
            Γ = 2 * ((V[3] - V[2]) / (S[3] - S[2]) - (V[2] - V[1]) / (S[2] - S[1])) / (S[3] - S[1])
        end
        if n == 0
            thetaval2 = V[1]
            Θ = (thetaval1 - thetaval2) / (4 * Δt)
            Δ_adjusted = Δ + (μ - rt + σ_bs^2 / 2) * Δt * s * Γ
            return (value, Δ, Γ, Θ, Δ_adjusted)
        end
    end
end

function crr_full_simple(optType::String, S0::Union{Float64,Int}, K::Union{Float64,Int}, r::Float64, σ::Float64, T::Union{Float64,Int}, N::Int)
    #=
    let's first implement and check the extended binomial tree Pelsser and Vorst greek calculations described in Cruz-Dias paper
    as per figure 3 in that paper, we need to extend our binomial tree:  not only do we need the S(0,0) root and it's descendants
    but we need to start the tree at S(-2,-1) = S/(U*D) as the root and then build the tree from there so we can approximate delta by
    Δ = (V(0,1)) - V(0,-1)) / (S(0,1) -  S(0,-1)) and
    Γ = 2[(V(0,1) - V(0,0) / (S(0,1) - S(0,0)) - (V(0,0) - V(0,-1)) / (S(0,0) - S(0,-1))] / (S(0,1) - S(0,-1))
    =#
    if lowercase(optType)[1] == 'c'
        mult = 1
    else
        mult = -1
    end
    Δt = T / (N - 2)  # we need to extend the tree by 2 steps to get the right delta, but we don't want to increase the stdev on the terminal slice
    #Δt = T / N
    R = exp(r * Δt)
    U = exp(σ * √Δt)
    D = 1 / U
    p = (R - D) / (U - D)
    q = (U - R) / (U - D)
    S = zeros(N + 1, N + 1)
    nr = N + 1
    nc = N + 1
    for r = nr:-1:1
        for c = 1:r
            S[r, c] = S0 / (U * D) * U^(c - 1) * D^(r - c)  # when r=3, c=2, this is S0 / (U*D) * U * D = S0
        end
    end
    V = 0 .* S
    V[nr, :] = max.(0, mult .* (S[nr, :] .- K))
    for r = nr-1:-1:1
        for c = 1:r
            #x = mult * (S[r, c] - K)
            y = (q * V[r+1, c] + p * V[r+1, c+1]) / R
            V[r, c] = max(x, y)
        end
    end
    return (V[3, 2], (V[3, 3] - V[3, 1]) / (S[3, 3] - S[3, 1]),
        2 * ((V[3, 3] - V[3, 2]) / (S[3, 3] - S[3, 2]) - (V[3, 2] - V[3, 1]) / (S[3, 2] - S[3, 1])) / (S[3, 3] - S[3, 1]))
    #return (V[1, 1], (V[3, 3] - V[3, 1]) / (S[3, 3] - S[3, 1]), 2 * ((V[3, 3] - V[3, 2]) / (S[3, 3] - S[3, 2]) - (V[3, 2] - V[3, 1]) / (S[3, 2] - S[3, 1])) / (S[3, 3] - S[3, 1]))
end

function cev_full_simple(optType::String, s::Union{Float64,Int}, K::Union{Float64,Int}, rt::Float64, q::Union{Float64,Int}, σ_bs::Float64, T::Union{Float64,Int}, N::Int, β::Float64)
    #=
    let's first implement and check the extended binomial tree Pelsser and Vorst greek calculations described in Cruz-Dias paper
    as per figure 3 in that paper, we need to extend our binomial tree:  not only do we need the S(0,0) root and it's descendants
    but we need to start the tree at S(-2,-1) = S/(U*D) as the root and then build the tree from there so we can approximate delta by
    Δ = (V(0,1)) - V(0,-1)) / (S(0,1) -  S(0,-1)) and
    Γ = 2[(V(0,1) - V(0,0) / (S(0,1) - S(0,0)) - (V(0,0) - V(0,-1)) / (S(0,0) - S(0,-1))] / (S(0,1) - S(0,-1))
    =#
    if lowercase(optType)[1] == 'c'
        mult = 1
    else
        mult = -1
    end
    Δt = T / (N - 2)  # we need to extend the tree by 2 steps to get the right delta, but we don't want to increase the stdev on the terminal slice
    R = exp(rt * Δt)
    nr = N + 1
    nc = N + 1
    #σ * S^β = σ_bs * S, so σ = σ_bs * S^(1-β)
    σ = σ_bs * s^(1 - β)
    α = 1 - β

    function x(S, σ, α)
        S^α / (α * σ)
    end

    function xInv(x, σ, α)
        # x = S^α / (α * σ)
        # so S = (x * α * σ)^(1/α)
        if x > 0
            return (x * α * σ)^(1 / α)
        else
            return 0
        end
    end

    function p_up(S, σ, α, rt, q, Δt)
        xval = x(S, σ, α)
        S_up = xInv(xval + √Δt, σ, α)
        S_down = xInv(xval - √Δt, σ, α)
        arg = (S * exp((rt - q) * Δt) - S_down) / (S_up - S_down)

        if xval > 0 && 0 <= arg <= 1
            return arg
        elseif xval <= 0 || arg < 0
            return 0
        elseif xval > 0 && arg > 1
            return 1
        end


    end

    X = zeros(nr, nc) .+ NaN
    S = zeros(nr, nc) .+ NaN
    S[1, 1] = s
    X[1, 1] = x(S[1, 1], σ, α)
    for r = 2:nr
        X[r, 1] = X[r-1, 1] - √Δt
        S[r, 1] = xInv(X[r, 1], σ, α)
        for c = 2:r
            X[r, c] = X[r, c-1] + 2√Δt
            S[r, c] = xInv(X[r, c], σ, α)
        end
    end

    V = 0 .* S
    V[nr, :] = max.(0, mult .* (S[nr, :] .- K))
    for r = nr-1:-1:1
        for c = 1:r
            p = p_up(S[r, c], σ, α, rt, q, Δt)
            contVal = ((1 - p) * V[r+1, c] + p * V[r+1, c+1]) / R
            exVal = mult * (S[r, c] - K)
            V[r, c] = max(exVal, contVal)
        end
    end
    return (V[3, 2], (V[3, 3] - V[3, 1]) / (S[3, 3] - S[3, 1]),
        2 * ((V[3, 3] - V[3, 2]) / (S[3, 3] - S[3, 2]) - (V[3, 2] - V[3, 1]) / (S[3, 2] - S[3, 1])) / (S[3, 3] - S[3, 1]))
    #or, per Hull's book: return (V[1, 1], (V[3, 3] - V[3, 1]) / (S[3, 3] - S[3, 1]), 2 * ((V[3, 3] - V[3, 2]) / (S[3, 3] - S[3, 2]) - (V[3, 2] - V[3, 1]) / (S[3, 2] - S[3, 1])) / (S[3, 3] - S[3, 1]))
end

#=
function conversion_am_simple(treeType, S, K, μ, r, q, σ, T, N)
    if lowercase(treeType) == "crr"
        fn = crr_efficient_simple
    elseif lowercase(treeType) == "trinomial"
        fn = trinomial_simple
    end

    callPrice = fn("c", S, K, μ, r, q, σ, T, N)[1]
    putPrice = fn("p", S, K, μ, r, q, σ, T, N)[1]
    conversion = putPrice - callPrice + S - K
    return (conversion)
end

function implied_repo_from_conversion_am_simple(convPrice, treeType, S, K, μ, r, σ, T, N)
    global max_allowable_q = r  #recall q = r - r_repo
    global min_allowable_q = r

    offset_vec = vcat(10, 0, -1:-1:-3, -10:-10:-100, -200:-100:-900, -1000:-1000:-3000) ./ 100
    pvec = zeros(length(offset_vec))
    qvec = r .- offset_vec
    for i = 1:length(qvec)
        pvec[i] = conversion_am_simple(treeType, S, K, μ, r, qvec[i], σ, T, N)
        if (i > 1) && ((pvec[i] - convPrice) * (pvec[i-1] - convPrice) < 0) # i.e. the soln lies in between rvec[i] and rvec[i-1]
            min_allowable_q = qvec[i-1]
            max_allowable_q = qvec[i]
            break
        elseif (pvec[i] > convPrice)
            return NaN # pvec is a decreasing fn of r so as we walk r in the negative direction, pvec > mktPrice means we can't cross from below and should exit
        end
    end
    if abs(pvec[end]) > 0
        return NaN # pvec[end] will be non-zero only if no soln was found
    end
    obj(x) = (conversion_am_simple(treeType, S, K, μ, r, x, σ, T, N) - convPrice)^2
    res = Optim.optimize(obj, min_allowable_q, max_allowable_q)
    if Optim.converged(res)
        return 100 * (r - Optim.minimizer(res))
    end
end
=#

#import ForwardDiff:Dual
#=
function lsmc_am(optType, S₀, K, r, σ, T, Nstep, Npath)

    Δt = T / Nstep
    Z = exp(-r * Δt)
    type = typeof(S₀ * exp((r-σ^2/2) * Δt + σ * √Δt * 0.1))
    S = Array{type}(undef, Nstep + 1, 2 * Npath)  #each path is a column of Nsteps
    
    Random.seed!(1234)
    rn = randn(Nstep,Npath)
    rn = hcat(rn, -rn)
    for p = 1:2*Npath
        S[1, p] = S₀
        for n = 2:Nstep + 1
            S[n, p] = S[n-1, p] * exp((r - σ^2 / 2) * Δt + σ * √Δt * rn[n-1,p])
        end
    end

    rn = nothing
    GC.gc()
    
    if lowercase(optType)[1] == "c"
        mult = 1.0
    else
        mult = -1.0
    end
    val = max.(mult .* (S .- K), 0.0)

    #following Longstaff and Schwartz, we'll use weighted Laguerre polynomials
    L₀(x) = 1
    L₁(x) = x
    L₂(x) = x^2 #1/2 * (3x^2 - 1)
    L₃(x) = x^3 #1/2 * (5x^3 - 3x)
    L₄(x) = x^4 #1/8 * (35x^4 - 30x^2 + 3)

    for t = Nstep:-1:2
        indx  = findall((@view val[t,:]) .!= 0) #which paths at step n are in the money
        Svec = @view S[t, indx]
        X = [L₀.(Svec) L₁.(Svec) L₂.(Svec) ]#L₃.(Svec) L₄.(Svec)]
        Y = zeros(size(X,1)) # y at step n is the discounted val of all steps > n (only 1 of which can be nonzero)
        rng = 1:(Nstep+1-t)
        df = Z.^collect(rng)
        Threads.@threads for i in  eachindex(Y)
            c = indx[i]
            Y[i] = sum(val[t+1:Nstep+1,c] .* df)
        end
        β = X \ Y
        #Ŷ = β[1] * X[:,1] + β[2] * X[:,2] + β[3] * X[:,3] + β[4] * X[:,4] + β[5] * X[:,5]
        Ŷ = sum(transpose(β .* X'), dims=2)
        Threads.@threads for i in eachindex(Y)
            c = indx[i]
            if val[t,c] > Ŷ[i]
                val[t+1:end,c] .= 0
                #val[r,c] has already been initialized to exercise value
            else
                val[t,c] = 0
            end
        end
    end

    df = Z.^collect(1:Nstep)
    runsum = 0
    for c in 1:size(val,2)
        runsum = runsum + sum(df .* val[2:end,c])
    end
    return max(mult .* (S₀ - K), runsum/size(val,2))
    
end
=#

function lsmc_am(optType::String, S₀::Union{Float64,Int}, K::Union{Float64,Int}, r::Union{Float64,Int}, σ::Float64, T::Union{Float64,Int}, Nstep::Int, Npath::Int)
    #using Loess
    Δt = T / Nstep
    Z = exp(-r * Δt)
    type = typeof(S₀ * exp((r - σ^2 / 2) * Δt + σ * √Δt * 0.1))

    #Random.seed!(1234)

    # let's work with s = S/K so the numbers are closer to 1
    if true
        s = Array{type}(undef, Nstep + 1, 2 * Npath)  #each path is a column of Nsteps
        s[1, :] .= S₀ / K
        rn = randn(Nstep, Npath)
        rn = hcat(rn, -rn)   
        for p = 1:2*Npath, n = 2:Nstep+1
            s[n, p] = s[n-1, p] * exp((r - σ^2 / 2) * Δt + σ * √Δt * rn[n-1, p])
        end
        rn = nothing
        GC.gc()
    else
        s = Array{type}(undef, Nstep + 1, Npath)  #each path is a column of Nsteps
        s[1, :] .= S₀ / K
        for p = 1:Npath, n = 2:Nstep+1
           s[n, p] = s[n-1, p] * exp((r - σ^2 / 2) * Δt + σ * √Δt * randn())
        end
    end 


    if lowercase(optType)[1] == 'c'
        mult = 1.0
    else
        mult = -1.0
    end
  
    #following Longstaff and Schwartz, we'll use weighted Laguerre polynomials
    L₀(x) = exp(-x / 2)
    L₁(x) = exp(-x / 2) * (-x + 1)
    L₂(x) = exp(-x / 2) * (x^2 - 4x + 2) /2
    L₃(x) = exp(-x / 2) * (-x^3 + 9x^2 - 18x + 6) / 6
    L₄(x) = exp(-x / 2) * (x^4 - 16x^3 + 72x^2 - 96x + 24) / 24

    exVal = max.(mult .* (s[end,:] .- 1.0), 0.0)
    contVal = copy(exVal)

    for t = Nstep:-1:2
        exVal = max.(mult .* (s[t, :] .- 1.0), 0.0)
        contVal = Z * contVal
        indx = findall(exVal .> 10 * eps()) #which paths at step n are in the money
        if length(indx) .> 0
            svec = @view s[t, indx]
            X = [L₀.(svec) L₁.(svec) L₂.(svec) L₃.(svec) L₄.(svec)]
            Y = contVal[indx] * Z
            β = X \ Y
            Ŷ = sum(transpose(β .* X'), dims=2)
            #model = loess(svec, Y, span = 0.5, degree = 4)
            #Ŷ = predict(model, svec)
            Threads.@threads for i in eachindex(indx)
                p = indx[i]
                if exVal[p] > Ŷ[i]
                    contVal[p] = exVal[p]
                end
            end
        end
    end

    return max(mult * (S₀ - K), Z * K * mean(contVal))

end


#import ForwardDiff:Dual
#lsmc_am_put(Dual(100, 1, 0), 90, 0.05, Dual(0.3, 0, 1), 180 / 365, 1000, 10000)
#will return (value, delta, vega)

# do corrected deltas and prices as per Wilmott pg 773
# add mu to inputs and add delta_adjusted to outputs
# note that setting the new input mu = r will make vol_adjusted = vol
# for reference:
# Δ_adjusted = Δ + (μ - r + σ^2/2) * Δt * S * Γ ... > Δ if μ > r
# σ_adjusted = σ * (1 + Δt * (μ - r) * (r - μ - σ^2) / (2 * σ^2))  ... < σ if μ > r

function cev_simple_with1Jump_stepwiseLocalVol(optType, S, K, rt, q, dt_trade, dt_expiry, N, σ=[0.25], dt_vol=[today()], dt_jump=today(), U=1, D=1, β=0.5)
    ##  THIS DOESN'T WORK YET - THE ISSUE IS THAT S ON FINAL SLICE NODES IS COMPUTED AS XINVERSE OF X ON FINAL NODES
    ##  AND X(N, i) IS COMPUTED AS X0 PLUS JUMPS THAT ONLY DEPEND ON A SINGLE STEP SIZE AND VOL
    ##  ***  WORK TO BE DONE!  ***
    if lowercase(optType)[1] == 'c'
        mult = 1
    else
        mult = -1
    end
    #Δt = T / (N - 2)  # we need to extend the tree by 2 steps to get the right delta, but we don't want to increase the stdev on the terminal slice
    #R = exp(rt * Δt)
    #σ * S^β = σ_bs * S, so σ = σ_bs * S^(1-β)
    #σ = σ_bs * S^(1 - β)
    α = 1 - β

    puj = (1 - D) / (U - D)
    # make sure the vol inputs go at least as far into the future as dt_expiry
    if dt_vol[end] < dt_expiry
        dt_vol[end] = dt_expiry
    end
    dates = sort(unique(vcat(dt_trade, dt_vol, dt_expiry)))
    dates = dates[dates.<=dt_expiry]
    # let's create a dataframe to hold the date-dependent info we need
    dum = DataFrame(startDate=dates[1:end-1], endDate=dates[2:end], σ=NaN .+ zeros(length(dates) - 1), var=NaN .+ zeros(length(dates) - 1), Δt=NaN .+ zeros(length(dates) - 1))
    for i = 1:nrow(dum)
        indx = findall(dt_vol .>= dum.endDate[i])
        dum.σ[i] = σ[indx[1]]
        dum.var[i] = (dum.endDate[i] - dum.startDate[i]).value / 365 * dum.σ[i]^2
    end
    N = round(Int64, 2 * sum(dum.var) / minimum(dum.var)) # at least 2 steps in each bucket
    N = max(N, (dt_expiry - dt_trade).value)
    #N = ceil(Int64, 365 * sum(dum.var)/minimum(dum.σ)^2) # at least 1step per day in all buckets
    varPerStep = sum(dum.var) / N
    dum.Δt = varPerStep ./ ((dum.σ) .^ 2)
    dum.u .= exp(√varPerStep)
    dum.d = 1 ./ dum.u
    dum.R = exp.((rt - q) * dum.Δt)
    dum.p_up = (dum.R - dum.d) ./ (dum.u - dum.d)
    dum.p_down = (dum.u - dum.R) ./ (dum.u - dum.d)
    dum.i_start = zeros(Int64, nrow(dum))
    dum.i_end = zeros(Int64, nrow(dum))
    dum.i_start[1] = 0
    dum.i_end[1] = round(Int64, (dum.endDate[1] - dum.startDate[1]).value / 365 / dum.Δt[1])
    for i = 2:nrow(dum)
        dum.i_start[i] = dum.i_end[i-1]
        dum.i_end[i] = dum.i_start[i] + round(Int64, (dum.endDate[i] - dum.startDate[i]).value / 365 / dum.Δt[i])
    end
    N = dum.i_end[end] #redefine N in case the rounding screwed things up
    # now we have to go back and fix time steps
    dum.Δt = Dates.value.(dum.endDate .- dum.startDate) ./ (dum.i_end .- max.(0, dum.i_start)) / 365

    tStepVec = []
    for i = 1:nrow(dum)
        nstep = dum.i_end[i] - max(0, dum.i_start[i])
        tStepVec = vcat(tStepVec, dum.Δt[i] * ones(nstep))
    end
    tVec = cumsum(tStepVec)
    tVec = vcat(0, tVec)

    dtVec = dum.startDate[1] .+ Day.(ceil.(365 * tVec))
    if dt_jump <= dt_expiry
        #############   THIS IS A KLUDGE FOR WHEN Δt SPACING CAUSES dt_jump TO BE SKIPPED OVER
        row = findall(dum.startDate .< dt_jump .&& dum.endDate .>= dt_jump)[1]
        i_jump = round(Int64, min(dum.i_end[row], dum.i_start[row] + ceil((dt_jump - dum.startDate[row]).value / 365 / dum.Δt[row])))
        #i_jump = findmin(abs.(Dates.value.(dtVec .- dt_jump)))[2]
    end

    function x(S, σ, α)
        S^α / (α * σ)
    end

    function xInv(x, σ, α)
        # x = S^α / (α * σ)
        # so S = (x * α * σ)^(1/α)
        if x > 0
            return (x * α * σ)^(1 / α)
        else
            return 0
        end
    end

    function p_up(S, σ, α, rt, q, Δt)
        xval = x(S, σ, α)
        S_up = xInv(xval + √Δt, σ, α)
        S_down = xInv(xval - √Δt, σ, α)
        arg = (S * exp((rt - q) * Δt) - S_down) / (S_up - S_down)

        if xval > 0 && 0 <= arg <= 1
            return arg
        elseif xval <= 0 || arg < 0
            return 0
        elseif xval > 0 && arg > 1
            return 1
        end
    end

    S0 = S
    x0 = x(S0, σ, α)
    # populate terminal values
    X = zeros(N + 1)
    S = zeros(N + 1)
    V = zeros(N + 1)

    if dt_jump > dt_expiry
        n = N
        for i = 0:n
            X[i+1] = x0 + sqrt(dum.Δt[end]) * (-n + 2 * i)
            S[i+1] = xInv(X[i+1], σ[end], α)
            V[i+1] = max(0, mult * (S[i+1] - K))
        end
        # propagate backward
        for n = N-1:-1:0
            r = findall(dum.i_start .<= n .&& dum.i_end .> n)[1]
            for i = 0:n
                X[i+1] = x0 + sqrt(dum.Δt[r]) * (-n + 2 * i)
                S[i+1] = xInv(X[i+1], σ[r], α)
                p = p_up(S[i+1], σ[r], α, rt, q, dum.Δt[r])
                contVal = ((1 - p) * V[i+1] + p * V[i+2]) / dum.R[r]
                exVal = mult * (S[i+1] - K)
                V[i+1] = max(exVal, contVal)
            end
        end
        return V[1]
    else
        # up tree
        V_up = zeros(N + 1)
        n = N
        for i = 0:n
            X[i+1] = x0 + sqrt(dum.Δt[end]) * (-n + 2 * i)
            S[i+1] = xInv(X[i+1], σ[end], α)
            V_up[i+1] = max(0, mult * (U * S[i+1] - K))
        end
        for n = N-1:-1:i_jump
            r = findall(dum.i_start .<= n .&& dum.i_end .> n)[1]
            for i = 0:n
                X[i+1] = x0 + sqrt(dum.Δt[r]) * (-n + 2 * i)
                S[i+1] = xInv(X[i+1], σ[r], α)
                p = p_up(S[i+1], σ[r], α, rt, q, dum.Δt[r])
                contVal = ((1 - p) * V_up[i+1] + p * V_up[i+2]) / dum.R[r]
                exVal = mult * (U * S[i+1] - K)
                V_up[i+1] = max(exVal, contVal)
            end
        end
        # down tree
        V_down = zeros(N + 1)
        n = N
        for i = 0:n
            X[i+1] = x0 + sqrt(dum.Δt[end]) * (-n + 2 * i)
            S[i+1] = xInv(X[i+1], σ[end], α)
            V_down[i+1] = max(0, mult * (D * S[i+1] - K))
        end
        for n = N-1:-1:i_jump
            r = findall(dum.i_start .<= n .&& dum.i_end .> n)[1]
            for i = 0:n
                X[i+1] = x0 + sqrt(dum.Δt[r]) * (-n + 2 * i)
                S[i+1] = xInv(X[i+1], σ[r], α)
                p = p_up(S[i+1], σ[r], α, rt, q, dum.Δt[r])
                contVal = ((1 - p) * V_down[i+1] + p * V_down[i+2]) / dum.R[r]
                exVal = mult * (D * S[i+1] - K)
                V_up[i+1] = max(exVal, contVal)
            end
        end
        n = i_jump - 1
        ###### NEXT LINE BUGS
        r = findall(dum.i_start .<= n .&& dum.i_end .> n)[1]
        ###### PREV LINE BUGS
        for i = 0:n
            X[i+1] = x0 + sqrt(dum.Δt[r]) * (-n + 2 * i)
            S[i+1] = xInv(X[i+1], σ[r], α)
            p = p_up(S[i+1], σ[r], α, rt, q, dum.Δt[r])
            contVal = (puj * ((1 - p) * V_up[i+1] + p * V_up[i+2]) + (1 - puj) * ((1 - p) * V_down[i+1] + p * V_down[i+2])) / dum.R[r]
            exVal = mult * (S[i+1] - K)
            V[i+1] = max(exVal, contVal)
        end
        for n = i_jump-2:-1:0
            r = findall(dum.i_start .<= n .&& dum.i_end .> n)[1]
            for i = 0:n
                X[i+1] = x0 + sqrt(dum.Δt[r]) * (-n + 2 * i)
                S[i+1] = xInv(X[i+1], σ[r], α)
                p = p_up(S[i+1], σ[r], α, rt, q, dum.Δt[r])
                contVal = ((1 - p) * V[i+1] + p * V[i+2]) / dum.R[r]
                exVal = mult * (S[i+1] - K)
                V_up[i+1] = max(exVal, contVal)
            end
        end
        return V[1]
    end
end


function crr_simple_with1Jump(optType, S, K, r, q, σ, T, N, t_jump, U, D)
    # puj*U + (1-puj)*D = 1 <===> puj = (1-D)/(U - D) where U = multiple for up jump, D = multiple for down move, and puj = prob. of up jump
    # we will need at least 3 strikes at a given expiry to fit vol, p and U
    if lowercase(optType)[1] == 'c'
        mult = 1
    else
        mult = -1
    end
    Δt = T / (N - 2)
    ##### jump params
    i_jump = ceil(Int, 2 + N * t_jump / T)
    puj = (1 - D) / (U - D)
    #####
    R = exp((r - q) * Δt)
    u = exp(σ * √Δt)
    d = 1 / u
    p_up = (R - d) / (u - d)
    p_down = (u - R) / (u - d)

    thetaval1 = thetaval2 = Θ = value = Δ = Γ = 0.0
    V = zeros(N + 1)
    #
    V_up = zeros(N + 1)
    n = N
    for i = 0:n
        V_up[i+1] = max(0, mult * (S * U * u^i * d^(N - i) - K))
    end
    for n = N-1:-1:i_jump
        for i = 0:n
            x = mult * (S * U * u^i * d^(n - i) - K)
            y = (p_down * V_up[i+1] + p_up * V_up[i+2]) / R
            V_up[i+1] = max(x, y)
        end
    end
    V_down = zeros(N + 1)
    n = N
    for i = 0:n
        V_down[i+1] = max(0, mult * (S * D * u^i * d^(N - i) - K))
    end
    for n = N-1:-1:i_jump
        for i = 0:n
            x = mult * (S * D * u^i * d^(n - i) - K)
            y = (p_down * V_down[i+1] + p_up * V_down[i+2]) / R
            V_down[i+1] = max(x, y)
        end
    end
    n = i_jump - 1
    for i = 0:n
        x = mult * (S * u^i * d^(n - i) - K)
        y = (puj * (p_down * V_up[i+1] + p_up * V_up[i+2]) + (1 - puj) * (p_down * V_down[i+1] + p_up * V_down[i+2])) / R
        V[i+1] = max(x, y)
    end
    for n = i_jump-2:-1:0
        for i = 0:n
            x = mult * (S * u^i * d^(n - i) - K)
            y = (p_down * V[i+1] + p_up * V[i+2]) / R
            V[i+1] = max(x, y)
        end
        if n == 4
            thetaval1 = V[3]
        end
        if n == 2
            value = V[2]
            Δ = (V[3] - V[1]) / (S * u^2 - S * d^2)
            Γ = 2 * ((V[3] - V[2]) / (S * u^2 - S * u * d) - (V[2] - V[1]) / (S * u * d - S * d^2)) / (S * u^2 - S * d^2)
        end
        if n == 0
            thetaval2 = V[1]
            Θ = (thetaval1 - thetaval2) / (4 * Δt)
            return (value, Δ, Γ, Θ)
        end
    end
end

function crr_simple_with1Jump_stepwiseLocalVol(optType, S, K, R, q, dt_trade, dt_expiry, N, σ=[0.25], dt_vol=[today()], dt_jump=today(), U=1, D=1)
    if lowercase(optType)[1] == 'c'
        mult = 1
    else
        mult = -1
    end
    ndays = (dt_expiry - dt_trade).value
    T = ndays / 365
    ## prob of upward binary i_jump
    puj = (1 - D) / (U - D)

    #=
    we will make sure the nodes join up by varying the time step size to ensure equal variance per step as per 
    Hull section 21.5:
    Given V = <σ²>ₜT and nstep = N we choose each Δt such that σ²(tᵢ) Δtᵢ = V/N, i.e. tᵢ₊₁ = tᵢ + (V/N) / σ²(tᵢ)
    and we take U = exp(√(V/N)) and D=1/U.  As with constant vol, p_up(tᵢ) = (Rᵢ - d)/(u-d) where Rᵢ=exp((r-q)*Δtᵢ)
    =#

    # make sure the vol inputs go at least as far into the future as dt_expiry
    if dt_vol[end] < dt_expiry
        dt_vol[end] = dt_expiry
    end
    #
    dates = sort(unique(vcat(dt_trade, dt_vol, dt_expiry)))
    dates = dates[dates.<=dt_expiry]
    # let's create a dataframe to hold the date-dependent info we need
    dum = DataFrame(startDate=dates[1:end-1], endDate=dates[2:end], σ=NaN .+ zeros(length(dates) - 1), var=NaN .+ zeros(length(dates) - 1), Δt=NaN .+ zeros(length(dates) - 1))
    for i = 1:nrow(dum)
        indx = findall(dt_vol .>= dum.endDate[i])
        dum.σ[i] = σ[indx[1]]
        dum.var[i] = (dum.endDate[i] - dum.startDate[i]).value / 365 * dum.σ[i]^2
    end
    ## varPerStep will be sum(dum.var) / N and therefore dum.Δt will be varPerStep ./ (dum.σ).^2
    ## we need to make sure that dum.Δt <= 1 day so we need to change N so that sum(dum.var) / N / (dum.σ)^2
    ## is <=1 day for all intervals, i.e. that 1/365 >= sum(dum.var) / N / min(dum.σ^2)
    ## or, finally, N >= 365 * sum(dum.var) / min(dum.σ^2) 
    ## and then let's make it a bit bigger just for safety:  N = ceil(Int, 1.25*365 * sum(dum.var)/minimum(dum.σ)^2)

    #=
    ******************************
    we need dum.i_end[1] to be >=1

    let's make sure every vol bucket has at least 2 steps. this will minimize compute time (maximize step size) given the
    equal var / step constraint
    i.e. find N s.t.
    varPerStep = sum(dum.var) / N
    ∀ᵢ nstep in bucket i = Tᵢ / Δtᵢ = Tᵢ / (varPerStep / σᵢ²) =  N * σᵢ² * Tᵢ / sum(dum.var) >= 2 , i.e. N * varᵢ / totvar >= 2
    or, N >= 2 * totvar / varᵢ for all i
    i.e.  N >= 2 * totvar / minimium(varᵢ)
    #*******************************
    =#
    N = round(Int64, 2 * sum(dum.var) / minimum(dum.var)) # at least 2 steps in each bucket
    N = max(N, (dt_expiry - dt_trade).value)
    #N = ceil(Int64, 365 * sum(dum.var)/minimum(dum.σ)^2) # at least 1step per day in all buckets
    varPerStep = sum(dum.var) / N
    dum.Δt = varPerStep ./ ((dum.σ) .^ 2)
    dum.u .= exp(√varPerStep)
    dum.d = 1 ./ dum.u
    dum.R = exp.((R - q) * dum.Δt)
    dum.p_up = (dum.R - dum.d) ./ (dum.u - dum.d)
    dum.p_down = (dum.u - dum.R) ./ (dum.u - dum.d)
    dum.i_start = zeros(Int64, nrow(dum))
    dum.i_end = zeros(Int64, nrow(dum))
    dum.i_start[1] = 0
    dum.i_end[1] = round(Int64, (dum.endDate[1] - dum.startDate[1]).value / 365 / dum.Δt[1])
    for i = 2:nrow(dum)
        dum.i_start[i] = dum.i_end[i-1]
        dum.i_end[i] = dum.i_start[i] + round(Int64, (dum.endDate[i] - dum.startDate[i]).value / 365 / dum.Δt[i])
    end
    N = dum.i_end[end] #redefine N in case the rounding screwed things up
    # now we have to go back and fix time steps
    dum.Δt = Dates.value.(dum.endDate .- dum.startDate) ./ (dum.i_end .- max.(0, dum.i_start)) / 365

    tStepVec = []
    for i = 1:nrow(dum)
        nstep = dum.i_end[i] - max(0, dum.i_start[i])
        tStepVec = vcat(tStepVec, dum.Δt[i] * ones(nstep))
    end
    tVec = cumsum(tStepVec)
    tVec = vcat(0, tVec)

    dtVec = dum.startDate[1] .+ Day.(ceil.(365 * tVec))
    if dt_jump <= dt_expiry
        #############   THIS IS A KLUDGE FOR WHEN Δt SPACING CAUSES dt_jump TO BE SKIPPED OVER
        row = findall(dum.startDate .< dt_jump .&& dum.endDate .>= dt_jump)[1]
        i_jump = round(Int64, min(dum.i_end[row], dum.i_start[row] + ceil((dt_jump - dum.startDate[row]).value / 365 / dum.Δt[row])))
        #i_jump = findmin(abs.(Dates.value.(dtVec .- dt_jump)))[2]
    end

    thetaval1 = thetaval2 = Θ = value = Δ = Γ = 0.0
    V = zeros(N + 1)
    #
    if dt_jump > dt_expiry
        n = N
        u = dum.u[end]
        d = dum.d[end]
        for i = 0:n
            V[i+1] = max(0, mult * (S * u^i * d^(N - i) - K))
        end
        for n = N-1:-1:0
            r = findall(dum.i_start .<= n .&& dum.i_end .>= n + 1)[1]
            R = dum.R[r]
            p_up = dum.p_up[r]
            p_down = dum.p_down[r]
            u = dum.u[r]
            d = dum.d[r]
            for i = 0:n
                x = mult * (S * u^i * d^(n - i) - K)
                y = (p_down * V[i+1] + p_up * V[i+2]) / R
                V[i+1] = max(x, y)
            end
            if n == 4
                thetaval1 = V[3]
            end
            if n == 2
                value = V[2]
                Δ = (V[3] - V[1]) / (S * u^2 - S * d^2)
                Γ = 2 * ((V[3] - V[2]) / (S * u^2 - S * u * d) - (V[2] - V[1]) / (S * u * d - S * d^2)) / (S * u^2 - S * d^2)
            end
            if n == 0
                thetaval2 = V[1]
                Θ = (thetaval1 - thetaval2) / (4 * dum.Δt[1])
                return (value, Δ, Γ, Θ)
            end
        end
        return V[1]
    else
        V_up = zeros(N + 1)
        n = N
        u = dum.u[end]
        d = dum.d[end]
        for i = 0:n
            V_up[i+1] = max(0, mult * (S * U * u^i * d^(N - i) - K))
        end
        for n = N-1:-1:i_jump
            r = findall(dum.i_start .<= n .&& dum.i_end .>= n + 1)[1]
            R = dum.R[r]
            p_up = dum.p_up[r]
            p_down = dum.p_down[r]
            u = dum.u[r]
            d = dum.d[r]
            for i = 0:n
                x = mult * (S * U * u^i * d^(n - i) - K)
                y = (p_down * V_up[i+1] + p_up * V_up[i+2]) / R
                V_up[i+1] = max(x, y)
            end
        end
        V_down = zeros(N + 1)
        n = N
        u = dum.u[end]
        d = dum.d[end]
        for i = 0:n
            V_down[i+1] = max(0, mult * (S * D * u^i * d^(N - i) - K))
        end
        for n = N-1:-1:i_jump
            r = findall(dum.i_start .<= n .&& dum.i_end .>= n + 1)[1]
            R = dum.R[r]
            p_up = dum.p_up[r]
            p_down = dum.p_down[r]
            u = dum.u[r]
            d = dum.d[r]
            for i = 0:n
                x = mult * (S * D * u^i * d^(n - i) - K)
                y = (p_down * V_down[i+1] + p_up * V_down[i+2]) / R
                V_down[i+1] = max(x, y)
            end
        end
        n = i_jump - 1
        ###### NEXT LINE BUGS
        r = findall(dum.i_start .<= n .&& dum.i_end .> n)[1]
        ###### PREV LINE BUGS
        R = dum.R[r]
        p_up = dum.p_up[r]
        p_down = dum.p_down[r]
        u = dum.u[r]
        d = dum.d[r]
        for i = 0:n
            x = mult * (S * u^i * d^(n - i) - K)
            y = (puj * (p_down * V_up[i+1] + p_up * V_up[i+2]) + (1 - puj) * (p_down * V_down[i+1] + p_up * V_down[i+2])) / R
            V[i+1] = max(x, y)
        end
        for n = i_jump-2:-1:0
            r = findall(dum.i_start .<= n .&& dum.i_end .>= n + 1)[1]
            R = dum.R[r]
            p_up = dum.p_up[r]
            p_down = dum.p_down[r]
            u = dum.u[r]
            d = dum.d[r]
            for i = 0:n
                x = mult * (S * u^i * d^(n - i) - K)
                y = (p_down * V[i+1] + p_up * V[i+2]) / R
                V[i+1] = max(x, y)
            end
        end
        return V[1]
    end
end

#==#
function implied_vol_crr_simple_with1Jump_stepwiseLocalVol(optPrice, optType, S, K, R, q, dt_trade, dt_expiry, N, σ=[0.25], dt_vol=[today()], dt_jump=today(), U=1, D=1)
    # this is to imply the local vol for a fixed strike, assuming earlier local vols are known, jump params are known if jump is earlier
    # σ is an array with 1 element per expiry that goes up to the previous expiry
    # dt_vol is an array of expiries going up to the expiry at the end of the local vol bucket being solved for

    function obj(x)
        vol = copy(σ)
        indx = findall(isnan.(vol))
        vol[indx] .= x
        (crr_simple_with1Jump_stepwiseLocalVol(optType, S, K, R, q, dt_trade, dt_expiry, N, vol, dt_vol, dt_jump, U, D)[1] - optPrice)^2
    end
    res = Optim.optimize(obj, 0.1, 2, Brent(); show_trace=false)
    vol = Optim.minimizer(res)
    #res = Optim.optimize(obj, [σ[maxiumum(findall(σ .> 0))]], show_trace = false)
    #vol = Optim.minimizer(res)[1]
    #obj(x) = (crr_simple_with1Jump_stepwiseLocalVol(optType, S, K, R, q, dt_trade, dt_expiry, N, vcat(σ[1:end-1], x), dt_vol, dt_jump, U, D)[1] - optPrice)^2
    #res = Optim.optimize(obj, [0.5], NelderMead();Optim.Options(show_trace = true, abs_tol = .005))
    #vol = Optim.minimizer(res)
    vol
end

function implied_jump_params_crr_simple_with1Jump_stepwiseLocalVol(data, S, R, q, dt_trade, dt_expiry, N, σ=[0.25], dt_vol=[today()], dt_jump=today())
    # data is a dataframe with colums for midPrice, strike, optType
    # this is to imply U,D,σ for date range within which jump occurs

    function obj(x)
        vol = copy(σ)
        indx = findall(isnan.(vol))
        vol[indx] .= x[1]^2
        nr = nrow(data)
        myvals = zeros(nr)
        #Threads.@threads for i in 1:nr
        for i in 1:nr
            #myvals[i] = crr_simple_with1Jump_stepwiseLocalVol(data.optType[i], S, data.strike[i], R, q, dt_trade, dt_expiry, N, vcat(σ[1:end-1], 0 + x[1]^2), dt_vol, dt_jump, (1 + x[2]^2), 1 / (1 + x[3]^2))[1]   end
            myvals[i] = crr_simple_with1Jump_stepwiseLocalVol(data.optType[i], S, data.strike[i], R, q, dt_trade, dt_expiry, N, vol, dt_vol, dt_jump, (1 + x[2]^2), 1 / (1 + x[3]^2))[1]
        end
        return sum((myvals .- data.splinePrice) .^ 2)
    end
    initial_x = [sqrt(0.25), sqrt(0.2), 1]
    res = Optim.optimize(obj, initial_x; show_trace=false)
    if Optim.converged(res)
        vol = 100 * (0 + Optim.minimizer(res)[1]^2)
        U = 1 + Optim.minimizer(res)[2]^2
        D = 1 / (1 + Optim.minimizer(res)[3]^2)
        println("vol = $(round(vol))")
        println("U = $U")
        println("D = $D")
        return (vol, U, D)
    else
        println("Didn't converge :<")
        return (missing, missing, missing)
    end
end
#==#

function optionChainSmoother(optType, kVecIn, pxVecIn, S, dt_trade, dt_expiry, rt, q, kVecOut=kVecIn)
    #= we follow section 3 of Fengler, Arbitrage-free smoothing of the implied volatility surface
    which in turn is based on Green, P. J. and Silverman, B. W. (1994). Nonparametric regression and generalized linear models
    ... both are in the .../Documents/Option Stuff folder
    =#
    T = (dt_expiry - dt_trade).value / 365
    #make sure kVec is in ascending order and that pxVec matches...
    indx = sortperm(kVecIn)
    kVecIn = kVecIn[indx]
    pxVecIn = pxVecIn[indx]

    function splineFit(u, y, λ)
        n = length(u) - 2
        g = Vector{Float64}(undef, n)  # length n
        γ = Vector{Float64}(undef, n - 2) # length n-2
        h = diff(u)[2:end-1] # length n-1
        # because of the weird column labeling in Q (starting with 2)
        # let's define a standard matrix M with rows 1:n and columns 1:n-1
        # and then discard column 1 to define Q...
        M = zeros(n, n - 1) #label rows 1:n and label columns 2:n-1
        for j = 2:n-1
            M[j, j] = -h[j-1]^-1 - h[j]^-1
            M[j-1, j] = h[j-1]^-1
            M[j+1, j] = h[j]^-1
        end
        Q = M[:, 2:end]
        # we'll do a similar construction for R 
        M2 = zeros(n - 1, n - 1)
        for j = 2:n-1
            M2[j, j] = 1 / 3 * (h[j-1] + h[j])
        end
        for j = 2:n-2
            M2[j, j+1] = M2[j+1, j] = h[j] / 6
        end
        R = M2[2:end, 2:end]

        #define...
        yvec = vcat(y[2:end-1], zeros(n - 2))
        xvec = vcat(g, γ)
        A = [Q; -R']  # [n x n-2; n-2 x n-2] so A is 2n-2 x n-2
        B = zeros(2n - 2, 2n - 2)
        B[1:n, 1:n] = Diagonal(ones(n))

        B[n+1:2n-2, n+1:2n-2] = λ * R

        # we want to minimize over xvec the quadratic form -yvec'*xvec + 1/2 xvec' * B * xvec
        # subject to A' * xvec = 0


        model = Model(Ipopt.Optimizer)
        set_silent(model)
        @variable(model, xx[1:2n-2])
        @objective(model, Min, xx' * B * xx - 2 * (yvec' * xx))
        @constraint(model, A' * xx .== 0) # guarantees that spline values, 1st and 2nd derivs are continuous at knots
        @constraint(model, xx[n+1:end] .>= 0) # guarantees γᵢ > 0 so pdf >=0 everywhere
        if lowercase(optType) == "call"
            @constraint(model, xx[n] >= 0) # price positivity for highest strike call
            @constraint(model, xx[2:n] .- xx[1:n-1] .<= 0) # higher strike calls never worth more


        #@constraint(model, (xx[2] - xx[1]) >= -exp(rt * T) * (u[3] - u[2])) #monotonicity            
        #@constraint(model, xx[1] >= exp(-q * T) * S - exp(-rt * T) * u[2])
        #@constraint(model, xx[1] <= exp(-q * T) * S)

        elseif lowercase(optType) == "put"
            @constraint(model, xx[1] >= 0) # price positivity for lowest strike put
            @constraint(model, xx[2:n] .- xx[1:n-1] .>= 0) # higher strike puts never worth less


            #@constraint(model, (xx[2] - xx[1]) >= exp(rt * T) * (u[3] - u[2])) #monotonicity            
            #@constraint(model, xx[1] >= -exp(-q * T) * S + exp(-rt * T) * u[2])
            #@constraint(model, xx[1] <= exp(-rt * T) * u[2])

        end

        #
        #


        optimize!(model)
        solution_summary(model)
        x = value.(xx)
        g = x[1:n]
        γ = [0; x[n+1:end]; 0]
        return (g, γ, h)
    end

    function spline(û, u, g, γ, h)
        function spline_in(û)
            # find uᵢ, uᵢ₊₁ s.t. û ∈ [uᵢ, uᵢ₊₁]
            x = vcat(g, γ)
            i_above = max(3, findfirst(u .>= û))
            i_below = i_above - 1
            uᵢ = u[i_below]
            uᵢ₊₁ = u[i_above]
            hᵢ = h[i_below-1]
            gᵢ = x[i_below-1]
            gᵢ₊₁ = x[i_below]
            # recall that len(x) = 2n-2 and x = [g,γ] where γ is of length n-2 because it ignores 0,n
            γᵢ = γ[i_below-1]
            γᵢ₊₁ = γ[i_below]
            val = ((û - uᵢ) * gᵢ₊₁ + (uᵢ₊₁ - û) * gᵢ) / hᵢ - 1 / 6 * (û - uᵢ) * (uᵢ₊₁ - û) * ((1 + (û - uᵢ) / hᵢ) * γᵢ₊₁ + (1 + (uᵢ₊₁ - û) / hᵢ) * γᵢ)
            return val
        end
        if û <= u[end-1] && û >= u[2]
            return spline_in(û)
        elseif û > u[end-1]
            uₙ = u[end-1]
            uₙ₋₁ = u[end-2]
            gₙ = spline_in(uₙ)
            gₙ₋₁ = spline_in(uₙ₋₁)
            γₙ₋₁ = γ[end-1]
            slope = (gₙ - gₙ₋₁) / (uₙ - uₙ₋₁) + 1 / 6 * (uₙ - uₙ₋₁) * γₙ₋₁
            val = gₙ + (û - uₙ) * slope
            #println("û is above the interpolation interval")
            return val
        elseif û < u[2]
            u₁ = u[2]
            u₂ = u[3]
            g₁ = spline_in(u₁)
            g₂ = spline_in(u₂)
            γ₂ = γ[2]
            slope = (g₂ - g₁) / (u₂ - u₁) - 1 / 6 * (u₂ - u₁) * γ₂
            val = g₁ + (û - u₁) * slope
            #println("û is below the interpolation interval")
            return val
        end
    end

    λ = 100
    (g, γ, h) = splineFit(kVecIn, pxVecIn, λ)
    return (kVecOut, spline.(kVecOut, Ref(kVecIn), Ref(g), Ref(γ), Ref(h)))

end

function channelRule(px::Vector{Float64}, len::Int64, trailingStopSigs::Float64, volLen::Int64)
    # first compute trailing max and min (we'll put the max of S[r-len:r-1]) on row r of trailingMax, e.g.)
    # also compute trailing vol with same convention
    nd = length(px)
    rtn = vcat(NaN, diff(px) ./ px[1:end-1])
    σ = zeros(nd) .+ NaN
    for t in volLen+2:nd
        σ[t] = nanstd(rtn[t-volLen:t-1])
    end
    trailingMin = zeros(nd) .+ NaN
    trailingMax = zeros(nd) .+ NaN
    for t in len+1:nd
        trailingMax[t] = maximum(px[t-len:t-1])
        trailingMin[t] = minimum(px[t-len:t-1])
    end
    # we'll need arrays for position (1 or -1), pnl by step
    pos = zeros(nd)
    pnl = zeros(nd)
    global tradePnls = []
    tradePnl = 0
    for t in len+1:nd
        if px[t] > trailingMax[t]
            pos[t] = 1 # long entry condition is met
        elseif px[t] < trailingMin[t]
            pos[t] = -1 # short entry condition is met
        else
            pos[t] = pos[t-1] # position is carried forward from previous step
        end
        pnl[t] = pos[t-1] * (px[t] - px[t-1]) #this is where you'd subtract a friction term
        if isnan(pnl[t])
            pnl[t] = 0
        end
        if pos[t] == pos[t-1]
            tradePnl = tradePnl + pnl[t]
        else
            if tradePnl != 0 #i.e. we've exited prev trade and perhaps entered new one
                global tradePnls = vcat(tradePnls, tradePnl)
            end
            tradePnl = 0
        end
        if tradePnl < -trailingStopSigs * σ[t] * px[t-1]# trailing stop loss
            pos[t] = 0
            push!(tradePnls, tradePnl)
            tradePnl = 0
        end
    end
    return pnl, pos, σ, tradePnl, trailingMax, trailingMin
end

function ewma(px::Vector{Float64}, len::Int64)
    α = 1 - exp(-1 / len)
    nd = length(px)
    avgpx = zeros(nd) .+ NaN
    for t in 1:nd
        if t < len
            avgpx[t] = nanmean(px[1:t])
        elseif t >= len
            avgpx[t] = α * px[t] + (1 - α) * avgpx[t-1]
        end
    end
    return avgpx
end

function ewma_better(x::AbstractVector, len::Int64)
    # see ...Dropbox/Notes/Exponential Moving Average.pdf
    α = exp(-1 / len) # formula is ~ x̄ᵢ = (1-α) xᵢ + α x̄ᵢ₋₁
    n = length(x)
    x̄ = zeros(n)
    x̄[1] = x[1]
    for t in 2:n
        x̄[t] = (1 - α) / (1 - α^t) * x[t] + (1 - α^(t - 1)) / (1 - α^t) * α * x̄[t-1]
    end
    return x̄
end

function macd(px::Vector{Float64}, longLen::Int64, shortLen::Int64, signalLen::Int64)
    nd = length(px)
    macd_line = ewma(px, shortLen) - ewma(px, longLen)
    signal_line = ewma(macd_line, signalLen)
    pos = zeros(nd)
    pnl = zeros(nd)
    for t in longLen:nd
        if macd_line[t] > signal_line[t]
            pos[t] = 1
        elseif macd_line[t] < signal_line[t]
            pos[t] = -1
        end
        pnl[t] = pos[t-1] * (px[t] - px[t-1])
        if isnan(pnl[t])
            pnl[t] = 0
        end
    end
    return pnl, pos, macd_line, signal_line
end


"""
sample usage:
```
using Distributions, Plots
obs = rand(Normal(0, 1), 1000);
F⁰ = Normal(0, 1)
myqqplot(obs, F⁰, "Normal QQ-plot ")

```
"""
function myqqplot(obs, F⁰, title)
    nobs = length(obs)
    sort!(obs)
    quantiles⁰ = [quantile(F⁰, i / nobs) for i in 1:nobs]
    # Note that only n-1 points may be plotted, as quantile(F⁰,1) may be inf
    plot(quantiles⁰, obs, seriestype=:scatter, xlabel="Theoretical Quantiles", ylabel="Sample Quantiles", title=title, label="")
    plot!(obs, obs, label="")
end

function crr_yc_repovec(yc, optType, S, K, r_csa_sprd_pct, repoSpreadDF, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    # here repoSpreadDF is a dataframe with colNames = startDate, endDate, repoSprdPct
    (ycinfo, futs, intervals) = yc
    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end
    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)
    dt_exp_settle = workday(dt_expiry, exercise_settlement)

    # put in error trap code in case ycinfo.date[end] < dt_exp_settle

    ndays = (dt_expiry - dt_trade).value
    ndays2 = (dt_exp_settle - dt_trade).value
    T = ndays / 360
    N = ceil(Int, stepsPerDay * ndays) + 2
    N2 = ceil(Int, stepsPerDay * ndays2) + 2
    Δt = T / (N - 2)  # N.B. already in units of 360 day year

    #ycinfo.rateCC = ycinfo.rate .+ NaN
    #for i = 1:nrow(ycinfo)-1
    #    ycinfo.rateCC[i] = 360 * 100 * log(1 + 0.01 * ycinfo.rate[i] / 360)
    #end
    rateCC = 360 * 100 * log.(1 .+ 0.01 * ycinfo.rate ./ 360)


    dateVec = dt_trade .+ Day.(ceil.(round.((-2:N2-2) * 360 * Δt, digits=1)))
    dateVecNext = vcat(dateVec[2:end], dateVec[end])
    i_stock_settle = maximum(findall(dateVec .== dt_stock_settle))
    i_option_settle = maximum(findall(dateVec .== dt_option_settle))
    baseRateVec = zeros(length(dateVec)) .+ NaN
    dates = unique(dateVec[dateVec.>=dt_trade])
    for d ∈ dates
        rng1 = findall(dateVec .== d)
        #baseRateVec[rng1] .= ycinfo.rateCC[findall(ycinfo.date .== d)[1]]
        baseRateVec[rng1] .= rateCC[findall(ycinfo.date .== d)[1]]
    end
    stIndx = findall(dateVec .== dt_trade)[1]
    csaVec = baseRateVec .+ r_csa_sprd_pct
    # so that we only discount back to the option settlement date
    csaVec[findall(dateVec .<= dt_option_settle)] .= 0
    #repoVec = baseRateVec .+ r_repo_sprd_pct
    repoSpreadDF = repoSpreadDF[setdiff(1:nrow(repoSpreadDF), findall(repoSpreadDF.startDate .>= dt_exp_settle)), :]
    repoInfoDF = DataFrame(date = dateVec, baseRate = baseRateVec, repoSprd = zeros(length(dateVec)))
    for i in 1:nrow(repoSpreadDF)
        rng = findall(repoInfoDF.date .>= repoSpreadDF.startDate[i] .&& repoInfoDF.date .<= repoSpreadDF.endDate[i])
        repoInfoDF.repoSprd[rng] .= repoSpreadDF.repoSprdPct[i]
    end
    repoVec = Vector(repoInfoDF.baseRate) .+ Vector(repoInfoDF.repoSprd)
    # so that repo accrual begins on stock_settle date
    repoVec[findall(dateVec .<= dt_stock_settle)] .= 0

    dfCsaVec = ones(length(dateVec))
    dfCsaVec[stIndx:end] = [1; exp.(-cumsum(0.01 * csaVec[stIndx:end-1] * Δt))]
    dfRepoVec = ones(length(dateVec))
    dfRepoVec[stIndx:end] = [1; exp.(-cumsum(0.01 * repoVec[stIndx:end-1] * Δt))]
    aVec = exp.(0.01 * repoVec * Δt)
    Rvec = exp.(0.01 * csaVec * Δt)
    # so that we only discount back to the option settlement date
    Rvec[findall(dateVecNext .<= dt_option_settle)] .= 1

    σ = vol_pct / 100
    divfvs = zeros(Float64, N + 1)
    if sum(divamt) > 0
        # we want at each slice divfv = value of all divs w/ exdiv not in the past as of that slice
        i_exdiv = zeros(Int, length(divamt))
        i_divpmt = zeros(Int, length(divamt))
        for i = 1:length(divamt)
            i_divpmt[i] = floor(Int, 4 + stepsPerDay * ((dt_divpmt[i] - dt_trade).value - 1))
            i_exdiv[i] = floor(Int, 4 + stepsPerDay * ((dt_exdiv[i] - dt_trade).value - 1))
            # this gets discounting right even when divpmt falls between slices because stepsPerDay < 1
            #divfvs[1:i_exdiv[i]-1] .+= (divamt[i] * dfRepoVec[i_divpmt[i]]) ./ dfRepoVec[1:i_exdiv[i]-1]
            divfvs[1:i_exdiv[i]-1] .+= (divamt[i] * df(dt_divpmt[i], yc)) ./ dfCsaVec[1:i_exdiv[i]-1]
        end
    end
    # we will discount using R_csa from i_option_settle to end
    # we will grow S0 at R_repo from i_stock_settle to end
    # this means we'll go to the LAST slice on dt_option_settle and dt_stock_settle
    #i_option_settle = 1 + stepsPerDay * (dt_option_settle - dt_trade).value
    #i_stock_settle = 1 + stepsPerDay * (dt_stock_settle - dt_trade).value
    S0 = S - divfvs[1]
    σ = σ * S / S0

    U = exp(σ * √Δt)
    D = 1.0 / U
    pVec = (aVec .- D) / (U - D)
    # so that we only start repo accrual after stock settle
    pVec[findall(dateVecNext .<= dt_stock_settle)] .= (1 - D) / (U - D)
    qVec = 1 .- pVec

    local Z_csa::Float64
    local Z_repo::Float64
    local V::Vector{Float64}
    local p2use::Float64
    local q2use::Float64
    local R2use::Float64
    local px::Float64
    local x::Float64
    local y::Float64
    local tau::Float64

    if lowercase(optType)[1] == 'p'
        dt = dt_expiry
        dt_settle = workday(dt, exercise_settlement)
        i_settle = maximum(findall(dateVec .== dt_settle))
        Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, N + 1)]
        Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, N + 1)]
        #Z_csa = dfCsaVec[i_settle] / dfCsaVec[N+1]
        #Z_repo = dfRepoVec[i_settle] / dfRepoVec[N+1]

        V = Z_csa * [max(0.0, K - S0 * U^i * D^(N - i) / Z_repo) for i = 0:N]
        for n = N-1:-1:0
            dt = dt_trade + Day(ceil((n - 2) / (N - 2) * ndays))
            dt_settle = workday(dt, exercise_settlement)
            i_settle = maximum(findall(dateVec .== dt_settle))
            Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, n + 1)]
            Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, n + 1)]
            #Z_csa = dfCsaVec[i_settle] / dfCsaVec[n+1]
            #Z_repo = dfRepoVec[i_settle] / dfRepoVec[n+1]

            if dt == bd(dt)
                for i = 0:n
                    y = (qVec[n+1] * V[i+1] + pVec[n+1] * V[i+2]) / Rvec[n+1]
                    px = S0 * U^i * D^(n - i) + divfvs[n+1]
                    x = Z_csa * (K - px / Z_repo)
                    V[i+1] = max(x, y)
                end
            else
                for i = 0:n
                    V[i+1] = (qVec[n+1] * V[i+1] + pVec[n+1] * V[i+2]) / Rvec[n+1]
                end
            end
            if n == 2
                delta = (V[3] - V[1]) / (S0 * U^2 - S0 * D^2)
                gamma = 2 * ((V[3] - V[2]) / (S0 * U^2 - S0 * U * D) - (V[2] - V[1]) / (S0 * U * D - S0 * D^2)) / (S0 * U^2 - S0 * D^2)
                value = V[2]
                return (value, delta, gamma)
            end
        end
    elseif lowercase(optType)[1] == 'c'
        dt = dt_expiry
        dt_settle = workday(dt, exercise_settlement)
        i_settle = maximum(findall(dateVec .== dt_settle))
        Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, N + 1)]
        Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, N + 1)]
        #Z_csa = dfCsaVec[i_settle] / dfCsaVec[N+1]
        #Z_repo = dfRepoVec[i_settle] / dfRepoVec[N+1]
        V = Z_csa * [max(0.0, S0 * U^i * D^(N - i) / Z_repo - K) for i = 0:N]
        for n = N-1:-1:0
            dt = dt_trade + Day(ceil((n - 2) / (N - 2) * ndays))
            dt_settle = workday(dt, exercise_settlement)
            i_settle = maximum(findall(dateVec .== dt_settle))
            Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, n + 1)]
            Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, n + 1)]
            #Z_csa = dfCsaVec[i_settle] / dfCsaVec[n+1]
            #Z_repo = dfRepoVec[i_settle] / dfRepoVec[n+1]
            if dt == bd(dt)
                for i = 0:n
                    y = (qVec[n+1] * V[i+1] + pVec[n+1] * V[i+2]) / Rvec[n+1]
                    px = S0 * U^i * D^(n - i) + divfvs[n+1]
                    x = Z_csa * (px / Z_repo - K)
                    V[i+1] = max(x, y)
                end
            else
                for i = 0:n
                    V[i+1] = (qVec[n+1] * V[i+1] + pVec[n+1] * V[i+2]) / Rvec[n+1]
                end
            end
            if n == 2
                delta = (V[3] - V[1]) / (S0 * U^2 - S0 * D^2)
                gamma = 2 * ((V[3] - V[2]) / (S0 * U^2 - S0 * U * D) - (V[2] - V[1]) / (S0 * U * D - S0 * D^2)) / (S0 * U^2 - S0 * D^2)
                value = V[2]
                return (value, delta, gamma)
            end
        end
    end
end

function trinomial_yc_repovec(yc, optType, S, K, r_csa_sprd_pct, repoSpreadDF, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    # here r_repo_sprd_pct is a dataframe with colNames = startDate, endDate, repoSprdPct
    (ycinfo, futs, intervals) = yc
    premium_settlement = 1
    if dt_trade < Date(2024, 5, 28)
        stock_settlement = exercise_settlement = 2
    else
        stock_settlement = exercise_settlement = 1
    end
    dt_option_settle = workday(dt_trade, premium_settlement)
    dt_stock_settle = workday(dt_trade, stock_settlement)
    dt_exp_settle = workday(dt_expiry, exercise_settlement)

    # put in error trap code in case ycinfo.date[end] < dt_exp_settle

    ndays = (dt_expiry - dt_trade).value
    ndays2 = (dt_exp_settle - dt_trade).value
    T = ndays / 360
    N = ceil(Int, stepsPerDay * ndays) + 1
    N2 = ceil(Int, stepsPerDay * ndays2) + 1
    Δt = T / (N - 1)  # N.B. already in units of 360 day year

    #ycinfo.rateCC = ycinfo.rate .+ NaN
    #for i = 1:nrow(ycinfo)-1
    #    ycinfo.rateCC[i] = 360 * 100 * log(1 + 0.01 * ycinfo.rate[i] / 360)
    #end
    rateCC = 360 * 100 * log.(1 .+ 0.01 * ycinfo.rate ./ 360)

    dateVec = dt_trade .+ Day.(ceil.(round.((-1:N2-1) * 360 * Δt, digits=1)))
    dateVecNext = vcat(dateVec[2:end], dateVec[end])
    i_stock_settle = maximum(findall(dateVec .== dt_stock_settle))
    i_option_settle = maximum(findall(dateVec .== dt_option_settle))
    baseRateVec = zeros(length(dateVec)) .+ NaN
    dates = unique(dateVec[dateVec.>=dt_trade])
    for d ∈ dates
        rng1 = findall(dateVec .== d)
        #baseRateVec[rng1] .= ycinfo.rateCC[findall(ycinfo.date .== d)[1]]
        baseRateVec[rng1] .= rateCC[findall(ycinfo.date .== d)[1]]
    end
    stIndx = findall(dateVec .== dt_trade)[1]
    csaVec = baseRateVec .+ r_csa_sprd_pct
    # so that we only discount back to the option settlement date
    csaVec[findall(dateVecNext .<= dt_option_settle)] .= 0
    
    #repoVec = baseRateVec .+ r_repo_sprd_pct
    repoSpreadDF = repoSpreadDF[setdiff(1:nrow(repoSpreadDF), findall(repoSpreadDF.startDate .>= dt_exp_settle)), :]
    repoInfoDF = DataFrame(date=dateVec, baseRate=baseRateVec, repoSprd=zeros(length(dateVec)))
    for i in 1:nrow(repoSpreadDF)
        rng = findall(repoInfoDF.date .>= repoSpreadDF.startDate[i] .&& repoInfoDF.date .<= repoSpreadDF.endDate[i])
        repoInfoDF.repoSprd[rng] .= repoSpreadDF.repoSprdPct[i]
    end
    repoVec = Vector(repoInfoDF.baseRate) + Vector(repoInfoDF.repoSprd)
    # so that repo accrual begins on stock_settle date
    repoVec[findall(dateVec .<= dt_stock_settle)] .= 0

    dfCsaVec = ones(length(dateVec))
    dfCsaVec[stIndx:end] = [1; exp.(-cumsum(0.01 * csaVec[stIndx:end-1] * Δt))]
    dfRepoVec = ones(length(dateVec))
    dfRepoVec[stIndx:end] = [1; exp.(-cumsum(0.01 * repoVec[stIndx:end-1] * Δt))]
    Rvec = exp.(0.01 * csaVec * Δt)

    σ = vol_pct / 100
    divfvs = zeros(Float64, N + 1)
    if sum(divamt) > 0
        # we want at each slice divfv = value of all divs w/ exdiv not in the past as of that slice
        i_exdiv = zeros(Int, length(divamt))
        i_divpmt = zeros(Int, length(divamt))
        for i = 1:length(divamt)
            i_divpmt[i] = floor(Int, 3 + stepsPerDay * ((dt_divpmt[i] - dt_trade).value - 1))
            i_exdiv[i] = floor(Int, 3 + stepsPerDay * ((dt_exdiv[i] - dt_trade).value - 1))
            # this gets discounting right even when divpmt falls between slices because stepsPerDay < 1
            #divfvs[1:i_exdiv[i]-1] .+= (divamt[i] * dfRepoVec[i_divpmt[i]]) ./ dfRepoVec[1:i_exdiv[i]-1]
            divfvs[1:i_exdiv[i]-1] .+= (divamt[i] * df(dt_divpmt[i], yc)) ./ dfCsaVec[1:i_exdiv[i]-1]
        end
    end
    S0 = S - divfvs[1]
    σ = σ * S / S0

    u = exp(σ * sqrt(3 * Δt))
    m = 1.0
    d = 1.0 / u
    p_u_vec = 1 / 6 .+ (0.01 * repoVec .- σ^2 / 2) * sqrt(Δt / (12 * σ^2))

    local V::Vector{Float64}
    local px::Float64
    local x::Float64
    local y::Float64
    local Z_csa::Float64
    local Z_repo::Float64
    local i_settle::Int32

    if lowercase(optType)[1] == 'p'
        dt = dt_expiry
        dt_settle = workday(dt, exercise_settlement)
        i_settle = maximum(findall(dateVec .== dt_settle))
        i_settle = findall(dateVec .== dt_settle)[end]
        Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, N + 1)]
        Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, N + 1)]
        V = Z_csa * [max(0.0, K - S0 * u^i / Z_repo) for i = -N:N]
        for n = N-1:-1:0
            dt = dt_trade + Day(ceil((n - 1) / (N - 1) * ndays))
            dt_settle = workday(dt, exercise_settlement)
            i_settle = maximum(findall(dateVec .== dt_settle))
            i_settle = findall(dateVec .== dt_settle)[end]
            Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, n + 1)]
            Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, n + 1)]
            if dt == bd(dt)
                for i = -n:n
                    y = ((1 - p_u_vec[n+1] - 2 / 3) * V[i+n+1] + 2 / 3 * V[i+n+2] + p_u_vec[n+1] * V[i+n+3]) / Rvec[n+1]
                    px = S0 * u^i + divfvs[n+1]
                    x = Z_csa * (K - px / Z_repo)
                    V[i+n+1] = max(x, y)
                end
            else
                for i = -n:n
                    V[i+n+1] = ((1 - p_u_vec[n+1] - 2 / 3) * V[i+n+1] + 2 / 3 * V[i+n+2] + p_u_vec[n+1] * V[i+n+3]) / Rvec[n+1]
                end
            end
            if n == 1
                delta = (V[3] - V[1]) / (S0 * u - S0 * d)
                gamma = 2 * ((V[3] - V[2]) / (S0 * u - S0 * m) - (V[2] - V[1]) / (S0 * m - S0 * d)) / (S0 * u - S0 * d)
                value = V[2]
                return (value, delta, gamma)
            end
        end
    elseif lowercase(optType)[1] == 'c'
        dt = dt_expiry
        dt_settle = workday(dt, exercise_settlement)
        i_settle = maximum(findall(dateVec .== dt_settle))
        i_settle = findall(dateVec .== dt_settle)[end]
        Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, N + 1)]
        Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, N + 1)]
        V = Z_csa * [max(0.0, S0 * u^i / Z_repo - K) for i = -N:N]
        for n = N-1:-1:0
            dt = dt_trade + Day(ceil((n - 1) / (N - 1) * ndays))
            dt_settle = workday(dt, exercise_settlement)
            i_settle = maximum(findall(dateVec .== dt_settle))
            i_settle = findall(dateVec .== dt_settle)[end]
            Z_csa = dfCsaVec[i_settle] / dfCsaVec[max(i_option_settle, n + 1)]
            Z_repo = dfRepoVec[i_settle] / dfRepoVec[max(i_stock_settle, n + 1)]
            if dt == bd(dt)
                for i = -n:n
                    y = ((1 - p_u_vec[n+1] - 2 / 3) * V[i+n+1] + 2 / 3 * V[i+n+2] + p_u_vec[n+1] * V[i+n+3]) / Rvec[n+1]
                    px = S0 * u^i + divfvs[n+1]
                    x = Z_csa * (px / Z_repo - K)
                    V[i+n+1] = max(x, y)
                end
            else
                for i = -n:n
                    V[i+n+1] = ((1 - p_u_vec[n+1] - 2 / 3) * V[i+n+1] + 2 / 3 * V[i+n+2] + p_u_vec[n+1] * V[i+n+3]) / Rvec[n+1]
                end
            end
            if n == 1
                delta = (V[3] - V[1]) / (S0 * u - S0 * d)
                gamma = 2 * ((V[3] - V[2]) / (S0 * u - S0 * m) - (V[2] - V[1]) / (S0 * m - S0 * d)) / (S0 * u - S0 * d)
                value = V[2]
                return (value, delta, gamma)
            end
        end
    end
end

function option_am_yc_repovec(yc, treeType, optType, S, K, r_csa_sprd_pct, repoSpreadDF, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    if lowercase(treeType) == "crr"
        fn = crr_yc_repovec
    elseif lowercase(treeType) == "trinomial"
        fn = trinomial_yc_repovec
    end
    return fn(yc, optType, S, K, r_csa_sprd_pct, repoSpreadDF, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
end

function conversion_am_yc_repovec(yc, treeType, S, K, r_csa_sprd_pct, repoSpreadDF, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    callInfo = option_am_yc_repovec(yc, treeType, "c", S, K, r_csa_sprd_pct, repoSpreadDF, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    putInfo = option_am_yc_repovec(yc, treeType, "p", S, K, r_csa_sprd_pct, repoSpreadDF, vol_pct, dt_trade, dt_expiry, divamt, dt_exdiv, dt_divpmt, stepsPerDay)
    value = putInfo[1] - callInfo[1] + S - K
    delta = putInfo[2] - callInfo[2] + 1
    gamma = putInfo[3] - callInfo[3]
    return (value, delta, gamma)
end

function implied_dividend_from_conversion_am_yc_repovec(convPrice, yc, treeType, S, K, r_csa_sprd_pct, repoSpreadDF, vol_pct, dt_trade, dt_expiry, dt_exdiv, dt_divpmt, stepsPerDay, maxDiv2Px)
    n = length(dt_exdiv)
    obj(x) = (conversion_am_yc_repovec(yc, treeType, S, K, r_csa_sprd_pct, repoSpreadDF, vol_pct, dt_trade, dt_expiry, x * ones(n), dt_exdiv[1:n], dt_divpmt[1:n], stepsPerDay)[1] - convPrice)^2
    ans = NaN
    maxrng = 0.01
    while (maxrng < maxDiv2Px)
        res = Optim.optimize(obj, S * (maxrng - 0.01), S * maxrng, abs_tol=0.001)
        if Optim.converged(res) && Optim.minimum(res) < 0.0001
            ans = Optim.minimizer(res)
            break
        else
            maxrng = maxrng + 0.01
        end
    end
    return ans
end

"""
    chatterjee_correlation(x::AbstractVector, y::AbstractVector)

Calculate Chatterjee's correlation coefficient between two vectors.
"""
function chatterjee_correlation(x::AbstractVector, y::AbstractVector)
    if length(x) != length(y)
        throw(ArgumentError("Input vectors must have equal length"))
    end

    n = length(x)
    pairs = collect(zip(x, y))
    sort!(pairs, by=first)
    sorted_y = last.(pairs)

    ranks = zeros(Int, n)
    for i in 1:n
        subsequence = @view sorted_y[i:end]
        ranks[i] = searchsortedfirst(sort(subsequence), sorted_y[i])
    end

    rank_diffs = abs.(ranks[2:end] .- ranks[1:end-1])
    sum_diffs = sum(rank_diffs)

    return 1 - (3 * sum_diffs) / (n * n - n)
end

# this is to display a dataframe nicely in pluto
function display_df_scrollable(df::DataFrame;
    title::String="my title",
    max_height::Int=400,
    max_rows::Union{Int,Nothing}=nothing,
    showRowLabels::Bool=true,
    scrollable::Bool=true)

    #define a formatter to add commas
    function ptf(v, i, j)
        if !(v isa String)
            return replace(string(Int(round(v))), r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
        else
            return v
        end
    end

    # If max_rows is specified, limit the displayed rows
    display_df = isnothing(max_rows) ? df : first(df, max_rows)

    # Get the basic table HTML
    if showRowLabels == true
        table_html = pretty_table(String, display_df,
            backend=Val(:html),
            formatters=ptf,
            title=title,
            title_alignment=:C,
            row_labels=1:nrow(display_df),
            show_subheader=false,
            standalone=false,
            tf=tf_html_default)
    else
        table_html = pretty_table(String, display_df,
            backend=Val(:html),
            formatters=ptf,
            title=title,
            title_alignment=:C,
            show_subheader=false,
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
        """<caption style = "text-align: center;">""" => """<caption style="text-align: center; font-family: Georgia, serif; font-size: 24px; font-weight: bold; color: #2c3e50; padding: 10px; $(caption_position) background: white; z-index: 1; ">""") #border:none;

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
