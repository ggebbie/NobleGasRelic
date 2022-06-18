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

# get age distribution
g = agedistribution.(loc)
tg = taudeltaresponse()

Δg = g[1] - g[2]

# PyPlot version, not currently showing
figure(2)
clf()
line1, = PyPlot.plot(tg,g[1],"black",label="35°N, 152°W, 3.5 km")
line2, = PyPlot.plot(tg,g[2],"red",label="20°S, 152°W, 3.5 km")
line3, = PyPlot.plot(tg,Δg,"green",label="Δ")
grid("true")
xlabel("Lag, τ [yr]")
ylabel("mass fraction per yr [1/yr]")
legend()
PyPlot.savefig(plotsdir("deltaresponse_NPACvSPAC.png"))

# try to use Plots 
# Plots.plot(tg,g[1],color=:black,label="35°N, 152°W, 3.5 km")
# Plots.plot!(tg,g[2],color=:red,label="20°S, 152°W, 3.5 km")
# Plots.plot!(tg,g[1]-g[2],color=:green,label="Δ")
# plot!(xlabel="Lag, τ [yr]",ylabel="mass fraction per yr [1/yr]")
# Plots.savefig(plotsdir("deltaresponse_NPACvSPAC.png"))

# underdetermined Gauss-Markov solution
#tmp =
tₚ = 2022 .- tg

# prior statistics of solution
σₓ = 5 *ones(length(tₚ))   # assume 20 mbar year-to-year variations

σcentury = 20;

# permit centennial-scale correlation in SLP
nt = length(tₚ)
Clong = Matrix{Float64}(undef,nt,nt)
T = 500 # timescale of correlation
for xx = 1:nt
    for yy = 1:nt
        Clong[xx,yy] = σcentury^2 * exp(-(tₚ[xx]-tₚ[yy])^2/T^2)
    end
end

Cₓₓ = Diagonal(σₓ.^2) + Clong

# Gauss-Markov solution method
Eᵀ = Δg
E = transpose(Eᵀ)
σy = 0.2 # assume ΔNe has error of 0.2 mbar (Jenkins pers.comm.)
x̃ = Cₓₓ*Eᵀ*((E*Cₓₓ*Eᵀ + σy^2)\ΔNe)

# reduction in uncertainty
P⁻ = Cₓₓ*Eᵀ*((E*Cₓₓ*Eᵀ + σy^2)\(E*Cₓₓ))
P  = Cₓₓ - P⁻

# what is uncertainty of DACP to LIA change
ΔDACPtoLIA = M*x̃
σDACPtoLIA = √(M*P*transpose(M))
σDACP = √(Mdacp*P*transpose(Mdacp))

# plot with reference to modern
imod = findall(x -> x ≥ 1979,tₚ)
ilia = findall(x -> tinterval[:LIA][1] ≤ x ≤ tinterval[:LIA][2],tₚ)
    
nt = length(tₚ)
Mmod = Matrix{Float64}(undef,nt,nt)
Mlia = Matrix{Float64}(undef,nt,nt)
# for each row, subtract modern value
for rr = 1:length(tₚ)
    Mmod[rr,rr] = 1
    Mmod[rr,imod] .-= 1 ./ length(imod)
    Mlia[rr,rr] = 1
    Mlia[rr,ilia] .-= 1 ./ length(ilia)
end


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

# max, min surface signal
tmax = tₚ[findmax(xₚ)[2]]
tmin = tₚ[findmin(xₚ)[2]]

# make diagnostic: size of Δ between max min, assuming timing doesn't change

imax = findall(x -> tmax - 100 ≤ x ≤ tmax + 100, tₚ)
imin = findall(x -> tmin - 100 ≤ x ≤ tmin + 100, tₚ)

# Magnitude diagnostic
M = zeros(1,length(tₚ))
M[imax] .= 1/length(imax)
M[imin] .= -1/length(imax)

Mdacp = zeros(1,length(tₚ))
Mdacp[imax] .= 1/length(imax)

mag = M*xₚ

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
