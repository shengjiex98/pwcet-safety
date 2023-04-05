
module Experiments

export SamplerPWCET

import Random

struct SamplerPWCET <: Random.Sampler{BitVector}
    p::Real
    H::Integer
end

function Random.rand!(rng::Random.AbstractRNG, a::BitVector, sp::SamplerPWCET)
    for i = 1:sp.H
        a[i] = Random.rand(rng) < sp.p
    end
    a
end

function Random.rand(rng::Random.AbstractRNG, sp::SamplerPWCET)
    a = falses(sp.H)
    Random.rand!(rng, a, sp)
end

end
