using Revise
using NobleGasRelic
using DrWatson
#using PyPlot
using LinearAlgebra
using TMI
using DataFrames
using Interpolations
#using PlotlyJS
using Plots
using OrderedCollections
using CSV

tinterval = OrderedDict(:MOD => (1860, 2022),
                :LIA => (1350,1860),
                :MCA => (800, 1350),
                :DACP => (400,800),
                #:DACP2 => (550,650),
                 :RWP => (-250,400),
                 :preRWP => (-Inf,-250))

longname = Dict(:MOD => "Modern Warming",
                :LIA => "Little Ice Age",
                :MCA => "Medieval Climate Anomaly",
                :DACP => "Dark Ages Cold Period",
                #:DACP2 => "Dark Ages Cold Period (strict)",
                :RWP => "Roman Warm Period",
                :preRWP => "Pre-Roman Warm Period")

# add dates to longname
longnamelabel = Dict{Symbol,String}()
for (kk,vv) in longname
    longnamelabel[kk] = longname[kk]*" "*string(tinterval[kk])*" CE"
end

vintage = collect(keys(tinterval))
depth = collect(2000:500:4000)

# lon must match grid exactly or else "dropped dims" error
lon = [162, 210]

params = @strdict vintage depth tinterval longnamelabel
dicts = dict_list(params)

# avoid error with NetCDF key already existing by
# moving existing file to SAVE version
isfile(datadir("vintages_TMI_4x4_2012.nc")) &&
    mv(datadir("vintages_TMI_4x4_2012.nc"),
       datadir("vintages_TMI_4x4_2012_SAVE.nc"),force=true)

# planviews
map(vintages_planview,dicts)

## sections
params = @strdict vintage lon tinterval longnamelabel
dicts = dict_list(params)
map(vintages_section,dicts)

# get 2 age distributions (TTDs), ultimately we have info about their difference
n = 2
loc = Vector{Tuple}(undef,n)
loc[1] = (360-152,35,3500) # North Pacific
loc[2] = (360-152,-10,3500) # South Pacific
compare_deltaresponses(loc)

# N.Pac: why more RWP than DACP water?
diagnose_deltaresponse(loc[1])

# try simple inversion
#    local g = vintagedistribution(tinterval[vintage][1],tinterval[vintage][2],Δ,τ)

# set up a DataFrame
# ax1 = "location"
# ax2 = "longitude [°E]"
# ax3 = "latitude [°N]"
# ax4 = "depth [m]"

col1 = "Vintage"
col2 = "Vintage Name"
col2b = "Years CE"
col3 = "Southern Region"
col4 = "Northern Region"
col5 = "SLP Anomaly [mbar]"
col6 = "SLP Error [mbar]"

defs = Dict(col1 => vintage,
            col2 => [longnamelabel[vv] for vv in vintage],
            col2b => [tinterval[vv] for vv in vintage])
df = DataFrame(defs)

gnorth = Dict{Symbol,Float64}()
gsouth = Dict{Symbol,Float64}()
for vv in vintage
    println(vv)
    gtmp = vintage_atloc(vv,loc)
    gnorth[vv] = gtmp[1]
    gsouth[vv] = gtmp[2]
end

insertcols!(df, col3 => [round(100gsouth[vv],digits=1) for vv in vintage])
insertcols!(df, col4 => [round(100gnorth[vv],digits=1) for vv in vintage])

E = Matrix(df)[:,4:5]
ΔE = transpose(E[:,2]-E[:,1])/100

# find time difference between vintages
t̄ = Dict{Symbol,Float64}()
for vv in vintage
    t̄[vv] =  (tinterval[vv][1] + tinterval[vv][2])/2
end
t̄[:preRWP] = -500.0  # instead of Inf

# make a covariance matrix that penalizes differences
# greater than 1 mbar/century

# make a covariance matrix
nv = length(vintage)
#D = Matrix{Float64}(undef,nv,nv)
S⁻ = zeros(Float64,nv,nv)
scentury = 2 # mbar/ century expected trend
for (mm,ii) in enumerate(vintage)
    for (nn,jj) in enumerate(vintage)
        if ii != jj
            Δt = abs(t̄[ii] - t̄[jj])
            δ = zeros(nv)
            δ[mm] = 1.0
            δ[nn] = -1.0
            global S⁻ += 1/(scentury*Δt/100)^2 * (δ * transpose(δ))
        end
    end

    # add constraint that MOD equals zero (Reference)
    if ii == :MOD
        S⁻[mm,mm] = 1/(0.01^2) # within 0.01
    end
    
end

Cₓₓ = inv(S⁻)
# Solve it.
ΔNe = 3.5 # mbar
σΔNe = 0.5 # mbar

xtmp,P = gaussmarkovsolution(transpose(ΔE),ΔNe,σΔNe,Cₓₓ)

#xtmp = (transpose(ΔE)*W⁻*ΔE + S⁻) \ (transpose(ΔE)*W⁻*y)
ỹ = (transpose(E)*xtmp)/100
ñ = ΔE*xtmp - ΔNe

x̃ = Dict{Symbol,Float64}()
σx̃ = Dict{Symbol,Float64}()
for (mm,ii) in enumerate(vintage)
    x̃[ii] = xtmp[mm]
    σx̃[ii] = √P[mm,mm]
end

insertcols!(df, col5 => [round(x̃[vv],digits=1) for vv in vintage])
insertcols!(df, col6 => [round(σx̃[vv],digits=1) for vv in vintage])
CSV.write(datadir("sixvintages.csv"),df)

################################
# try with PlotlyJS instead
# doesn't work because legend outside plot axes
plotlyjs(markershape=:auto)
#plotlyjs()
#PlotlyJS.purge!(plt)
#clf()
Plots.plot()
layout = Layout(legend=attr(x=1800.0, y=10.0))
for vv in vintage
    println(t̄[vv])
    println(x̃[vv])
    plt = Plots.plot!((t̄[vv],x̃[vv]),yerr=σx̃[vv],label=longnamelabel[vv])
end
#Plots.plot!(xlabel="Calendar Year [CE]",ylabel="Δ(SLP) [mbar]",legend=:bottomleft)
#Plots.plot!(xlabel="Calendar Year [CE]",ylabel="Δ(SLP) [mbar]",)
#Plots.plot!(xlabel="Calendar Year [CE]",ylabel="Δ(SLP) [mbar]",Layout(legend=attr(yanchor="top", xanchor="right")))
Plots.plot!(xlabel="Calendar Year [CE]",ylabel="Δ(SLP) [mbar]",layout)
plotname = "sixvintages"
Plots.savefig(plotsdir(plotname*".svg"))

################################
# try with GR
# doesn't work because legend outside plot axes
gr()
#clf()
Plots.plot()
for vv in vintage
    println(t̄[vv])
    println(x̃[vv])
    plt = Plots.plot!((t̄[vv],x̃[vv]),yerr=σx̃[vv],label=longnamelabel[vv])
end
Plots.plot!(xlabel="Calendar Year [CE]",ylabel="Δ(SLP) [mbar]",legend=:bottomleft)
plotname = "sixvintages"
Plots.savefig(plotsdir(plotname*".svg"))
Plots.savefig(plotsdir(plotname*".pdf"))

### RIBBON PLOTS
gr()
plotlyjs()
Plots.plot()
t̄sort = zeros(nv)
x̃sort = zeros(nv)
σx̃sort = zeros(nv)
for (ii,vv) in enumerate(vintage)
    t̄sort[ii] = t̄[vv]
    x̃sort[ii] = x̃[vv]
    σx̃sort[ii] = σx̃[vv]
end
Plots.plot!(t̄sort,x̃sort,ribbon=σx̃sort,color=:grey,linealpha=0,legend=false,markershape=:none)
Plots.plot!([t̄sort'; t̄sort'], [(x̃sort-σx̃sort)';(x̃sort+σx̃sort)'],color=:black,linewidth=1)
Plots.plot!(xlabel="Calendar Year [CE]",ylabel="Δ(SLP) [mbar]",legend=false)
plotname = "sixvintages_ribbon"
Plots.savefig(plotsdir(plotname*".svg"))
Plots.savefig(plotsdir(plotname*".pdf"))

## add to the spreadsheet to put all results on there.
nr = nrow(df)
for vv in vintage
    colname1 = String(vv)*" [%]"
    colname2 = String(vv)*" SLP [mbar]"
    colname3 = String(vv)*" σ(SLP) [mbar]"
    #colname4 = String(vv)*" influence [mbar]"

    SLPvals = repeat([x̃[vv]],nr)
    insertcols!(df, colname1, colname2 => SLPvals, after=true)

    σvals= repeat([σx̃[vv]],nr)
    insertcols!(df, colname2, colname3 => σvals, after=true)

    #σvals= repeat([x̃[vv]],nr)
    #insertcols!(df, colname2, colname3 => σvals, after=true)

end
lastcolname1 = String(:preRWP)*" σ(SLP) [mbar]"
lastcolname2 = "Total SLP influence [mbar]"
insertcols!(df, lastcolname1, lastcolname2 => ỹ, after=true)

#dfmat = Matrix(df)
#dfdiff = Vector{Any}(undef,ncol(df))

#dfdiff[5:end] = dfmat[1,5:end] .- dfmat[2,5:end]
#dfdiff[1:4] .= ""
#insert!(df,
# insert 3rd row
#insert!.(eachcol(df), 3*ones(nc), dfdiff)
#df2 = DataFrame(dfdiff,names(df))



#Thomas R. Knutson, Jeffrey Ploshay. Journal of Climate. DOI: 10.1175/JCLI-D-19-0997.1: 2 mbar/century trends in CMIP5
#https://en.wikipedia.org/wiki/North_Atlantic_oscillation#/media/File:Winter-NAO-Index.svg: Annual change about 2 mbar/yr
# Wunsch 1999 \Psi(s) =  0.66s−0.22 mb 2/cycle/yr,
# power law = 0.2, wunsch 1990

# what is 1000 year expected variation?
Φ1000 = 0.66*(1/1000)^(-0.22)
Φ1 = 0.66*1^(-0.22)
A1000 = √(Φ1000/2)
A1 = √(Φ1/2)

# underdetermined Gauss-Markov solution
tₚ = 2022 .- tg

σcentury = 10
Tlong = 500

# prior statistics of solution
σₓ = 1  # assume 5 mbar year-to-year variations

Cₓₓ = priorcovariance(tₚ,σₓ,σcentury,Tlong)

ΔNe = 3.5 # mbar
σΔNe = 0.2 # assume ΔNe has error of 0.2 mbar (Jenkins pers.comm.)

x̃,P = gaussmarkovsolution(Δg,ΔNe,σΔNe,Cₓₓ)

F = trendmatrix(tₚ[end],tₚ[begin],tₚ)
trend,σtrend = NobleGasRelic.propagate(F,x̃,P)

# try various time intervals to compute trends
Tlist = 1800:200:2000

for T in Tlist

    ntrend = findall(x -> x == tₚ[end] + T, tₚ)
    ntrend = ntrend[1]
    trend = zeros(ntrend)
    σtrend = zeros(ntrend)
    ttrend = zeros(ntrend)

    for tt = 1:ntrend
        println(tt)
        F = trendmatrix(tₚ[end+1-tt],tₚ[end+1-tt-T],tₚ)
        trend[tt],σtrend[tt] = propagate(F,x̃,P)
        ttrend[tt] = (tₚ[end+1-tt] + tₚ[end+1-tt-T])/2
    end

    # try with PlotlyJS instead
    plotlyjs()
    Plots.plot(ttrend,T*trend,color=:black)
    Plots.plot!(ttrend,T*trend,ribbon=T*σtrend,color=:grey)
    Plots.plot!(xlabel="Calendar Year [CE]",ylabel=string(T)*"-yr Δ(SLP) [mbar]",legend=false)
    plotname = "SLPtrends_"*string(T)*"yr"
    Plots.savefig(plotsdir(plotname*".svg"))
    Plots.savefig(plotsdir(plotname*".pdf"))
end

Mmod = averagematrix(1860,2022,tₚ)
Mlia = averagematrix(1860,2022,tₚ)

# what is uncertainty of DACP to LIA change
ΔDACPtoLIA = M*x̃
σDACPtoLIA = √(M*P*transpose(M))
σDACP = √(Mdacp*P*transpose(Mdacp))

# plot with reference to modern
#imod = findall(x -> x ≥ 1979,tₚ)
#ilia = findall(x -> tinterval[:LIA][1] ≤ x ≤ tinterval[:LIA][2],tₚ)
# nt = length(tₚ)
# Mmod = Matrix{Float64}(undef,nt,nt)
# Mlia = Matrix{Float64}(undef,nt,nt)
# # for each row, subtract modern value
# for rr = 1:length(tₚ)
#     Mmod[rr,rr] = 1
#     Mmod[rr,imod] .-= 1 ./ length(imod)
#     Mlia[rr,rr] = 1
#     Mlia[rr,ilia] .-= 1 ./ length(ilia)
# end


# get error bars
Pmod = Mmod*(P*transpose(Mmod))
σPmod = .√(diag(Pmod))

Plia = Mlia*(P*transpose(Mlia))
σPlia = .√(diag(Plia))

# try with Plots instead
Plots.plot(tₚ,Mlia*x̃,color=:black,label="Gauß-Markov surface signal")
Plots.plot!(tₚ,Mlia*x̃ .+ σPlia,color=:grey,label="+1σ")
Plots.plot!(tₚ,Mlia*x̃ .- σPlia,color=:grey,label="-1σ")
plot!(xlabel="Calendar Year [CE]",ylabel="SLP - SLP(LIA)  [mbar]",legend=:bottomleft)
Plots.savefig(plotsdir("surfacesignal_gaussmarkov_errors.png"))

# try with PyPlot 
# Plots.plot(tₚ,Mlia*x̃,color=:green,label="Gauß-Markov surface signal")
# plot!(xlabel="Calendar Year [CE]",ylabel="SLP - SLP(LIA)  [mbar]",legend=:bottomleft)
# PyPlot.figure(130)
#plot(tg1,θ1,"black",label="35°N, 152°W, 3.5 km")
#plot(tg2,θ2,"red",label="20°S, 152°W, 3.5 km")
PyPlot.plot(tₚ,Mlia*x̃,"green",label="surface signal")
figure()
PyPlot.errorbar(tₚ,Mlia*x̃,yerr=σPlia)
grid("true")
xlabel("Lag, τ [yr]")
ylabel("signal []")
legend()
savefig(plotsdir("minimalsurfacesignal.png"))

# minimum-energy surface timeseries

#G*θ = 1
#J = (G*θ) ̇ (G*θ
F = svd(transpose(Δg),full=true)
ΔNe = 3.5 # mbar

# J
#θ1 = dot(g1,θ̃)
#θ2 = dot(g2,θ̃)


# particular solution 
xₚ = F.Vt[1,:]*((F.U'*ΔNe)/F.S[1])

PyPlot.figure(130)
#plot(tg1,θ1,"black",label="35°N, 152°W, 3.5 km")
#plot(tg2,θ2,"red",label="20°S, 152°W, 3.5 km")
PyPlot.plot(tₚ,xₚ,"green",label="minimum-energy surface signal")
grid("true")
xlabel("Lag, τ [yr]")
ylabel("signal []")
legend()
savefig(plotsdir("minimalsurfacesignal.png"))

# try with Plots instead
Plots.plot(tₚ,xₚ,color=:green,label="minimum-energy surface signal")
plot!(xlabel="Calendar Year [CE]",ylabel="sea level pressure [mbar]",legend=:bottomleft)
Plots.savefig(plotsdir("minimalsurfacesignal.png"))

# make diagnostic: size of Δ between max min, assuming timing doesn't change
# Magnitude diagnostic
# max, min surface signal
#tmax = tₚ[findmax(xₚ)[2]]
#tmin = tₚ[findmin(xₚ)[2]]
# make diagnostic: size of Δ between max min, assuming timing doesn't change
# Magnitude diagnostic
#imax = findall(x -> tmax - 100 ≤ x ≤ tmax + 100, tₚ)
#imin = findall(x -> tmin - 100 ≤ x ≤ tmin + 100, tₚ)
# M = zeros(1,length(tₚ))
# M[imax] .= 1/length(imax)
# M[imin] .= -1/length(imax)

# Mdacp = zeros(1,length(tₚ))
# Mdacp[imax] .= 1/length(imax)

mag = magnitude(xₚ,tₚ)

# what if all of imax has null space projection of 20?
σ = 20
x = xₚ
for nn in imax
    global x += σ .* F.Vt[nn+1,:]
end

# make LIA more negative
σ = 20
x = xₚ
for nn in imin
    global x -= σ .* F.Vt[nn+1,:]
end


# one random realization is enough for plotting
# assume that annual variations can be 20 mbar
σ = 20
x = xₚ
for nn = 1:length(tₚ)-1
    global x += σ * rand().*F.Vt[nn+1,:]
end
#end

# minimal surface signal with other solutions
Plots.plot(tₚ,x .- Mdacp*x,color=:grey,label="higher-energy surface signal")
Plots.plot!(tₚ,xₚ .- Mdacp*xₚ,color=:red,label="minimum-energy surface signal")
plot!(xlabel="Calendar Year [CE]",ylabel="sea level pressure relative to DACP [mbar]",legend=:bottomleft)
Plots.savefig(plotsdir("surfacesignals.png"))

# get a range of possible DACP to LIA variations
# nreal = 10000
# mag = Vector{Float64}(undef,nreal)
# for nn = 1:nreal
#     σ = 20
#     x = xₚ
#     for nn = 1:length(tₚ)-1
#         x += σ * rand().*F.Vt[nn+1,:]
#     end
#     mag[nn] = (M*x)[1]
# end
