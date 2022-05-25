using Revise, TMItransient, TMI, PacificNobleGasRelic, DrWatson

Δ,τ = read_stepresponse()

g = vintage(1850,2022,Δ,τ)

# gMOD = vintage(1860,2022, Δ,τ) # modern warming
# gLIA = vintage(1450,1860, Δ,τ) # little ice age
# gMCA = vintage(950, 1250, Δ,τ) # medieval climate anomaly
# gLALIA = vintage(550, 950, Δ,τ) # late antique little ice age
# gLALIA = vintage(550, 650, Δ,τ) # late antique little ice age (strict)
# gRWP = vintage(1,550, Δ,τ) # roman warm period

vintages = Dict(:MOD => (1860, 2022),
                :LIA => (1450,1850),
                :MCA => (950, 1250),
                :LALIA => (550,950),
                :LALIA2 => (550,650),
                :RWP => (1,550))

vlongnames = Dict(:MOD => "Modern Warming",
                :LIA => "Little Ice Age",
                :MCA => "Medieval Climate Anomaly",
                :LALIA => "Late Antique Little Ice Age",
                :LALIA2 => "Late Antique Little Ice Age (strict)",
                :RWP => "Roman Warm Period")

vintagelist = keys(vintages)
depth = 3000

for v in vintagelist
    println(v)
    tlabel = "Vintage: "* vlongnames[v] * " " * string(vintages[v]) * " CE, depth="*string(depth)*"m"
    fname = "vintage_"*string(v)*"_"*string(depth)*"m.png"

    local g = vintage(vintages[v][1],vintages[v][2],Δ,τ)

    # Plan view plots
    PacificNobleGasRelic.planviewplotcartopy(100g, 3000, 0:5:100, titlelabel=tlabel)
    mv(plotsdir("vintage.png"),plotsdir(fname),force=true)
    
end

# get some TTDs so that we can take difference of TTDs
