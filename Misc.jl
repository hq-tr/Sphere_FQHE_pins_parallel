module MiscRoutine

# Basis format conversions
bin2dex(config::BitVector) = findall(config) .- 1

dex2bin(state::Vector{Int}, N_orb::Int) = BitVector([i-1 in state for i in 1:N_orb]) 

function dec2bin(num::Integer, N_orb::Int)
	rawbin = BitVector(digits(num,base=2))
	for i in 1:(N_orb-length(rawbin))
		push!(rawbin,0)
	end
	return rawbin
end

function dec2binreverse(num::Integer, N_orb::Int)
	rawbin = BitVector(digits(num,base=2))
	for i in 1:(N_orb-length(rawbin))
		push!(rawbin,0)
	end
	return reverse(rawbin)
end

function dec2bin(num::String,N_orb::Int)
	return dec2bin(parse(Int,num), N_orb)
end

# Normalization coefficient on the sphere
sphere_coef(S,m) =   sqfactorial(S-m)/sqfactorial(S+m+1, 2S+1)

# Miscellaneous functions for LLL physics

findLz(root::BitVector) = sum(root .* collect(0:1:length(root)-1))

findLzsphere(root::BitVector, S::Float64) = sum(root .* collect(-S:1:S))

function findLzsphere(root::BitVector)
	S = (length(BitVector)-1.)/2.
	return findLzsphere(root)
end

sqfactorial(N) = prod(map(sqrt, 1:N)) 
# square root of N!, avoid overflow for large N and more efficient than sqrt(factorial(big(N)))
# (overflow starts at N=21)

sqfactorial(n,N) = prod(map(sqrt, n:N))


# Miscellaneous function for the torus
get_k_vector(m::Int,Nx::Int,Ny::Int) = (m÷Nx)/Nx,(m%Nx)/Ny # These are actually coefficients [m₁,m₂] such that k = m₁b₁+m₂b₂

export bin2dex, sqfactorial, dex2bin, findLz, findLzsphere, sphere_coef, dec2bin, dec2binreverse,get_k_vector
end