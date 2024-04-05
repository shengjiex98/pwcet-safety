using Serialization
using Distributions: Pareto, quantile

push!(LOAD_PATH, "../src")
using Experiments

PATH = "../data/nmc-dist"
OUTPUT_PATH = "../data-proxy"
OUTPUT_FILE = "nmc-dist"
BATCHSIZE = 1_000_000
H = 100 * 0.02
Q_VALUES = 0.01:0.01:0.99
DIST = Pareto(1.5, 0.01)

p = 0.99
i99 = round(Int64, BATCHSIZE * p)
i_low, i_high = find_intervals(BATCHSIZE, p, 0.05, centered=true)[1]

function get_quantile(q, h)
    filename = generate_filename(BATCHSIZE, q, h, th=16)
    if !isfile("$PATH/$filename.jls")
        @info "Parameters without valid data" b q h
        return Inf
    end
    data = deserialize("$PATH/$filename.jls")
    [data[i][2] for i in (i_low, i99, i_high)]
end

quantiles = map(Q_VALUES) do q
    h = quantile(DIST, q)
    get_quantile(q, h)
end |> stack

@info size(quantiles)

serialize("$OUTPUT_PATH/$OUTPUT_FILE.jls", quantiles)

full = vcat(reshape(Q_VALUES, 1, :), reshape(periods, 1, :), quantiles)'
writedlm("$OUTPUT_PATH/$OUTPUT_FILE.csv", full, ',')

