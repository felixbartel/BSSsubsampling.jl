using LinearAlgebra
include("../src/BSSsubsampling.jl")

m = 50
M = 200
b = 1.5


# equispaced Fourier matrix ##################################################
Y = [ exp(2im*π*k*x)/sqrt(M) for x in (0:M-1)/M, k in -m/2:m/2-1 ]

println("\n===== Testing the BSS algorithm =====")
idcs, s = BSSsubsampling.bss(Y, b; A = 1, B = 1)

λ = svdvals(Y'*Y)
println("initial bounds:")
println("A = ", minimum(λ), ", B = ", maximum(λ))

λ = svdvals(M/length(idcs)*Y[idcs, :]'*Y[idcs, :])
println("subsampled bounds:")
println("A = ", minimum(λ), ", B = ", maximum(λ))


# random Fourier matrix ##################################################
Y = [ exp(2im*π*k*x)/sqrt(M) for x in rand(M), k in -m/2:m/2-1 ]

println("\n===== Testing the BSSperp algorithm =====")
idcs, s = BSSsubsampling.bss_perp(Y, b)

λ = svdvals(Y'*Y)
println("initial bounds:")
println("A = ", minimum(λ), ", B = ", maximum(λ))

λ = svdvals(M/length(idcs)*Y[idcs, :]'*Y[idcs, :])
println("subsampled bounds:")
println("A = ", minimum(λ), ", B = ", maximum(λ))


# random cosine matrix ##################################################
Y = [ √2*cos(π*k*x)/sqrt(M) for x in rand(M), k in 0:m-1 ]

println("\n===== Testing the plainBSS algorithm =====")
idcs, s = BSSsubsampling.bss_plain(Y, b)

λ = svdvals(Y'*Y)
println("initial bounds:")
println("A = ", minimum(λ), ", B = ", maximum(λ))

λ = svdvals(M/length(idcs)*Y[idcs, :]'*Y[idcs, :])
println("subsampled bounds:")
println("A = ", minimum(λ), ", B = ", maximum(λ))
