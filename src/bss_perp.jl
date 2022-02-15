"""
    idcs, s = bssperp(Y, b)
    idcs, s = bssperp(Y, b; Δ = 1e-14, verbose = false)
  
# Input
 - `Y::Matrix{<:Number}`: frame
 - `b::Number`: target oversampling factor
 - `Δ::Number`: stability parameter

# Output
 - `idcs::Vector{Int}`: BSSperp-subsampled frame elements
 - `s::Vector{Float64}`: corresponding weights

"""
function bss_perp(Y::Matrix{<:Number}, b::Number; Δ::Number = 1e-14, verbose = false)::Tuple{Vector{Int}, Vector{Float64}}
    # orthogonalize
    Y = Matrix(qr(Y).Q)
    # apply BSS
    idcs, s = bss(Y, b; A = 1, B = 1, Δ = Δ, verbose = verbose)
    # return
    return idcs, s
end
