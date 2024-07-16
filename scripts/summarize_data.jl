using Serialization
using Distributions: Pareto, quantile
using DelimitedFiles

push!(LOAD_PATH, "../src")
print(Base.load_path())
using Experiments

# const SYS = "EWB"
# const DIST = Pareto(1.5, 0.001)
# const PATH = "../data/nmc-dist/$SYS/$DIST"
# const OUTPUT_FILE = "nmc-dist-$SYS-$DIST"

const JOB_ID = 41703162
const file_num = 1:30

const P = 0.99

for i in file_num
    PATH = "../data/mpc/$JOB_ID/$i"
    OUTPUT_FILE = "mpc-$JOB_ID-$i"

    OUTPUT_PATH = "../data-proxy"

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
