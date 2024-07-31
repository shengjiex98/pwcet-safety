using Serialization
using Distributions: Pareto, quantile
using DelimitedFiles

push!(LOAD_PATH, "$(@__DIR__)/../src")
using Experiments

# const SYS = "EWB"
# const DIST = Pareto(1.5, 0.001)
# const PATH = "../data/nmc-dist/$SYS/$DIST"
# const OUTPUT_FILE = "nmc-dist-$SYS-$DIST"

const JOB_ID = 44669211
const FILE_NUM = 1:100

const P = 0.99

for i in FILE_NUM
    # PATH = "$(@__DIR__)/../data/mpc/$JOB_ID/$i"
    PATH = "$(@__DIR__)/../data/mpc-grid/$JOB_ID/$i"
    OUTPUT_PATH = "$(@__DIR__)/../data-proxy/mpc-grid/y/"
    mkpath(OUTPUT_PATH)
    OUTPUT_FILE = "$JOB_ID-$i"

    function parsefile(filename::String)
        b_s, q_s, h_s = split(filename, "-")[1:3]
        b = Int64(parse(Float64, b_s[2:end]))
        q = parse(Float64, q_s[2:end])
        h = parse(Float64, h_s[2:end])
        return b, q, h
    end

    @info "Reading files in $PATH..."
    flush(stderr)

    FILES = readdir(PATH)
    BATCHSIZE =parsefile(FILES[1])[1]

    @info "$(length(FILES)) files found. Batch size is $BATCHSIZE"
    flush(stderr)

    I99 = round(Int64, BATCHSIZE * P)
    I_LOW, I_HIGH = find_intervals(BATCHSIZE, P, 0.05, centered=true)[1]

    proxy = map(FILES) do filename
        b, q, h = parsefile(filename)
        data = deserialize("$PATH/$filename")
        vcat([q, h], [data[i][2] for i in [I_LOW, I99, I_HIGH]])
    end |> stack |> transpose

    @info size(proxy)
    flush(stderr)
    # @info proxy

    @info "Writing results to files..."
    flush(stderr)
    
    serialize("$OUTPUT_PATH/$OUTPUT_FILE.jls", proxy)
    writedlm("$OUTPUT_PATH/$OUTPUT_FILE.csv", proxy, ',')
end
