using Revise, PacificNobleGasRelic, DrWatson

tinterval = Dict(:MOD => (1860, 2022),
                :LIA => (1350,1850),
                :MCA => (800, 1350),
                :DACP => (400,800),
                :DACP2 => (550,650),
                :RWP => (1,550))

longname = Dict(:MOD => "Modern Warming",
                :LIA => "Little Ice Age",
                :MCA => "Medieval Climate Anomaly",
                :DACP => "Dark Ages Cold Period",
                :DACP2 => "Dark Ages Cold Period (strict)",
                :RWP => "Roman Warm Period")

vintage = collect(keys(tinterval))
depth = collect(2000:500:4000)
lon = [180, 220]

params = @strdict vintage depth tinterval longname
dicts = dict_list(params)

# planviews
map(vintages_planview,dicts)

## sections

params = @strdict vintage lon tinterval longname
dicts = dict_list(params)


# get some TTDs so that we can take difference of TTDs
