"""
    idcs, s = bss_plain(Y, b)
    idcs, s = bss_plain(Y, b; A = 1, B = 1, Δ = 1e-14, verbose = false)
  
# Input
 - `Y::Matrix{<:Number}`: frame
 - `b::Number`: target oversampling factor
 - `Δ::Number`: stability parameter

# Output
 - `idcs::Vector{Int}`: BSS-subsampled frame elements
 - `s::Vector{Float64}`: corresponding weights

"""
function bss_plain(Y::Matrix{<:Number}, b::Number; Δ::Number = 1e-14, verbose = false)::Tuple{Vector{Int}, Vector{Float64}}
    M, m = size(Y) # initial frame size

    m < 10 && error("BSS: dimension of frame elements has to be at least 3.")
    b <= 1+10/m && error("BSS: b has to be at least 1+10/m.")

    # extend with ceil(αm) elements
    α = (b-1)/(1+2*b)
    # oversampling factor for bss
    b2 = 1/(1-2α)

    # create duplicate frame elements such that M/ceil(m*α) is an integer
    # in the paper this was done with zero vectors
    tmp = Int(ceil(m*α))
    if mod(M, tmp) != 0
        Y = vcat(Y, Y[1:tmp-mod(M, tmp), :])
        M_new = M +tmp-mod(M, tmp)
    else
        M_new = M
    end

    # construct extension
    Y = hcat(sqrt(tmp/M_new) * kron(I(tmp), ones(Int(M_new/tmp))), Y)
    # orthogonalize
    Y = Matrix(qr(Y).Q)

    # apply BSS
    idcs, s = bss(Y, b2; Δ = Δ, verbose = verbose)

    # remove the artificial duplicate frame elements
    idcs = mod.(idcs, M)
    idcs[idcs .== 0] .= M
    p = sortperm(idcs)
    idcs = idcs[p]
    s = s[p]
    idx = 1
    while idx != length(idcs)
        if idcs[idx] == idcs[idx+1]
            s[idx] += s[idx+1]
            deleteat!(idcs, idx+1)
            deleteat!(s, idx+1)
        else
            idx += 1
        end
    end
    
    return idcs, s
end
