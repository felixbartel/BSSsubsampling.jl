"""
    idcs, s = bss(Y, b)
    idcs, s = bss(Y, b; A = 1, B = 1, Δ = 1e-14, verbose = false)
  
# Input
 - `Y::Matrix{<:Number}`: frame
 - `b::Number`: target oversampling factor
 - `A::Number`: initial lower frame bound
 - `B::Number`: initial upper frame bound
 - `Δ::Number`: stability parameter

# Output
 - `idcs::Vector{Int}`: BSS-subsampled frame elements
 - `s::Vector{Float64}`: corresponding weights

"""
function bss(Y::Matrix{<:Number}, b::Number; A::Number = NaN, B::Number = NaN, Δ::Number = 1e-14, verbose = false)::Tuple{Vector{Int}, Vector{Float64}}
    M, m = size(Y) # initial frame size
    n = Int(ceil(b*m)) # target frame size
    if verbose
        println("BSS: dimension of frame elements:      m = ", m)
        println("BSS: initial number of frame elements: M = ", M)
        println("BSS: target number of me elements:     n = ", n)
    end

    if isnan(A) || isnan(B) # compute the frame bounds if they are not given
        λ = svdvals(Y'*Y)
        A, B = λ[end], λ[1]
    end
    
    κ = (B/A+1)/2 + sqrt((B/A+1)^2/4-1)
    if b < κ^2 
        error("BSS: Oversampling factor b should be bigger than κ² = ", κ^2, ".")
    end

    l = -m*sqrt(b)*κ/(1+Δ)             # lower spectral barrier
    u = m*(sqrt(b)+b)/(sqrt(b)-1)*B/A  # upper spectral barrier
    δ_l = 1/(1+Δ)                      # lower spectral barrier shift
    δ_u = (sqrt(b)+1)/(sqrt(b)-1)*B/A  # upper spectral barrier shift
    δ_l0 = δ_l
    δ_u0 = δ_u
    ε_l0, ε_u0 = 0, 0 # will be initialized in the loop

    idcs = []
    s = []
    Ak = zeros(typeof(Y[1,1]), m, m)
    for k in ( verbose ? ProgressBar(1:n) : 1:n )
        # 1. compute eigenvalues
        λ, V_Ak = eigen(Hermitian(Ak))
        if minimum(λ) < l; error("BSS: spectral barrier λmin < l not valid."); end
        if maximum(λ) > u; error("BSS: spectral barrier λmax > u not valid."); end
        # 2. compute lower and upper potentials
        ε_l = sum(1 ./ (λ .- l))
        ε_u = sum(1 ./ (u .- λ))
        if k == 1 # save the first iteration
            ε_l0 = ε_l
            ε_u0 = ε_u
        end
        # 3. put δ ...
        δl = 1 / (1/δ_l0 - κ*ε_l0 + κ*ε_l)
        δu = 1 / (1/δ_u0 + ε_u0 - ε_u)
        # 4. increment l and u
        l += δ_l
        u += δ_u
        # 5. compute the factors
        f_l = sum(1 ./ (λ .- l))
        f_u = sum(1 ./ (u .- λ))
        # 6. find the next index
        tmpl = 1 ./ (λ .- l)
        tmpu = 1 ./ (λ .- u)
        L, U = 0, 0
        i = 1
        p = randperm(M)
        #p = 1:M
        while i <= M
            Uy = abs.( V_Ak'*Y[p[i],:] ).^2
            L = sum(tmpl.*tmpl.*Uy) / (f_l-ε_l) - sum(tmpl.*Uy) |> real
            U = sum(tmpu.*tmpu.*Uy) / (ε_u-f_u) - sum(tmpu.*Uy) |> real
            #L = Y[i, :]' * (Ak-u*I)^-2 * Y[i, :] / (f_l-ε_l) - Y[i, :]' * (Ak-u*I)^-1 * Y[i, :] |> real
            #U = Y[i, :]' * (Ak-l*I)^-2 * Y[i, :] / (ε_u-f_u) - Y[i, :]' * (Ak-l*I)^-1 * Y[i, :] |> real
            L-U >= Δ/(2*m)*(1-1/sqrt(b)) && break
            i += 1
        end
        i > M && error("BSS: did not find a frame element fulfilling the update condition.")
        i = p[i]
        # 7. compute weight and new matrix
        t = 2/(L+U)
        Ak += t*Y[i,:]*Y[i,:]'
        # 8. update
        idx = findfirst( idcs .== i)
        if isnothing(idx) # new element
            append!(idcs, i)
            append!(s, t)
        else # we already have this element but update the weight
            s[idx] += t
        end
    end
    γ = (sqrt(b)+1)^2/((sqrt(b)-1)*(sqrt(b)-κ))
    if verbose
        println("BSS: lower frame bound guarantee: A = ", A)
        println("BSS: upper frame bound guarantee: B = ", B*γ*(1+Δ))
    end
    s *= (A/l+B*γ*(1+Δ)/u)/2

    return idcs, s
end
