using Revise, PacificNobleGasRelic, DrWatson

tinterval = Dict(:MOD => (1860, 2022),
                :LIA => (1450,1850),
                :MCA => (950, 1250),
                :LALIA => (550,950),
                :LALIA2 => (550,650),
                :RWP => (1,550))

longname = Dict(:MOD => "Modern Warming",
                :LIA => "Little Ice Age",
                :MCA => "Medieval Climate Anomaly",
                :LALIA => "Late Antique Little Ice Age",
                :LALIA2 => "Late Antique Little Ice Age (strict)",
                :RWP => "Roman Warm Period")

vintage = collect(keys(tinterval))
depth = collect(2000:500:4000)

params = @strdict vintage depth tinterval longname
dicts = dict_list(params)

map(vintage_diagnostics,dicts)

# get some TTDs so that we can take difference of TTDs
