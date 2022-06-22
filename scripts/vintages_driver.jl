using Revise, NobleGasRelic, DrWatson, PyPlot, LinearAlgebra

tinterval = Dict(:MOD => (1860, 2022),
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
for (kk,vv) in longname
    longname[kk] *= (" "*string(tinterval[kk])*" CE")
end

vintage = collect(keys(tinterval))
depth = collect(2000:500:4000)

# lon must match grid exactly or else "dropped dims" error
lon = [162, 210]

params = @strdict vintage depth tinterval longname
dicts = dict_list(params)

# planviews
map(vintages_planview,dicts)

## sections

params = @strdict vintage lon tinterval longname
dicts = dict_list(params)

map(vintages_section,dicts)

# get 2 age distributions (TTDs), ultimately we have info about their difference
n = 2
loc = Vector{Tuple}(undef,n)
loc[1] = (360-152,35,3500) # North Pacific
loc[2] = (360-152,-20,3500) # South Pacific
compare_deltaresponses(loc)



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
