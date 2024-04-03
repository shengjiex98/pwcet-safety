using Serialization

push!(LOAD_PATH, "../lib")
using Experiments

function get_quantile(path, b, q, h, p; cap=Inf)
    filename = generate_filename(b, q, h, th=16)
    if !isfile("$path/$filename.jls")
        @info "Parameters without valid data" b q h
        return Inf
    end
    min(deserialize("$path/$filename.jls")[round(Int64, b * p)][2], cap)
end

path = "../data/nmc-samenom"
b_large = 1_000_000
qs = append!(collect(0.1:0.05:0.9), [0.99, 0.999])
hs = collect(0.005:0.0025:0.06)
p = 0.99
