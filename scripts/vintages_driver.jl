using Revise, NobleGasRelic, DrWatson

tinterval = Dict(:MOD => (1860, 2022),
                :LIA => (1350,1860),
                :MCA => (800, 1350),
                :DACP => (400,800),
                #:DACP2 => (550,650),
                :RWP => (-250,400))

longname = Dict(:MOD => "Modern Warming",
                :LIA => "Little Ice Age",
                :MCA => "Medieval Climate Anomaly",
                :DACP => "Dark Ages Cold Period",
                #:DACP2 => "Dark Ages Cold Period (strict)",
                :RWP => "Roman Warm Period")

# add dates to longname
for (kk,vv) in longname
    longname[kk] *= (" "*string(tinterval[kk])*" CE")
end


vintage = collect(keys(tinterval))
depth = collect(2000:500:4000)
lon = [162, 208]

params = @strdict vintage depth tinterval longname
dicts = dict_list(params)

# planviews
map(vintages_planview,dicts)

## sections

params = @strdict vintage lon tinterval longname
dicts = dict_list(params)

map(vintages_section,dicts)

# get some TTDs so that we can take difference of TTDs
#locs = Vector{Tuple}(undef,2)
loc1 = (360-152,35,3500)
loc2 = (360-152,-20,3500)

wis1= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,1)
[wis1[i] = interpindex(loc1,Δ[1].γ) for i in 1]

wis2= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,1)
[wis2[i] = interpindex(loc2,Δ[1].γ) for i in 1]

Δloc1 = Vector{Float64}(undef,length(Δ))
for tt in 1:length(Δ)
    Δloc1[tt] = observe(Δ[tt],wis1,Δ[tt].γ)[1]
end

Δloc2 = Vector{Float64}(undef,length(Δ))
for tt in 1:length(Δ)
    Δloc2[tt] = observe(Δ[tt],wis2,Δ[tt].γ)[1]
end


tg1,g1 = PacificNobleGasRelic.deltaresponse(Δloc1,τ)
tg2,g2 = PacificNobleGasRelic.deltaresponse(Δloc2,τ)

figure(2)
clf()
line1, = plot(tg1,g1,"black",label="35°N, 152°W, 3.5 km")
line2, = plot(tg2,g2,"red",label="20°S, 152°W, 3.5 km")
line3, = plot(tg1,g1-g2,"green",label="Δ")
grid("true")
xlabel("Lag, τ [yr]")
ylabel("mass fraction per yr [1/yr]")
legend()

# smallest climate signal to give a value of 1.
#G*θ = 1
G = g1 - g2

#J = (G*θ) ̇ (G*θ

tmp = (transpose(G)*G)\1
θ̃ = G*tmp
# J

θ1 = dot(g1,θ̃)
θ2 = dot(g2,θ̃)

figure(13)
#plot(tg1,θ1,"black",label="35°N, 152°W, 3.5 km")
#plot(tg2,θ2,"red",label="20°S, 152°W, 3.5 km")
plot(θ̃,"green",label="minimal surface signal to produce: NPAC-SPAC=1")
grid("true")
xlabel("Lag, τ [yr]")
ylabel("signal []")
legend()
savefig(plotsdir("minimalsurfacesignal.png"))
