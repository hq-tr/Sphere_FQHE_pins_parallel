module HilbertSpaceGenerator
using Combinatorics
const epsilon = 1e-8
function fullhilbertspace(N_el::Integer, N_orb::Integer)::Vector{UInt64}
    States = UInt64[]
    indVec = collect(0:N_orb-1)
    occLists = combinations(indVec, N_el)
    for occList in occLists
        State = 0
        for i in occList
            State += 2^i
        end
        push!(States, State)
    end 
    sort!(States)
    return States
end

function fullhilbertspace(N_el::Integer, N_orb::Integer, L_z::Number)::Vector{UInt64}
    States = UInt64[]
    S = 0.5(N_orb-1)
    indVec = collect(0:N_orb-1)
    occLists = combinations(indVec, N_el)
    for occList in occLists
        State = 0
        for i in occList
            State += 2^i
        end
        if abs(L_z- findLZ(State, S))<epsilon
            push!(States, State)
        end
    end 
    sort!(States)
    return States
end

function findLZ(State::UInt64,S::Real)
    Lz = 0
    while State != 0
        lowest_one = State & -State        # get the lowest "1" bit
        position = floor(Int, log2(lowest_one))
        Lz += -position+S
        State &= State - 1                # clear the lowest "1" bit
    end
    return Lz
end


function get_Lz_sphere(basis,coef,n_orb)
    S = (n_orb-1)/2
    Lz_list = map(x->findLZ(x,S), basis)
    return sum(Lz_list .* abs2.(coef))
end

export fullhilbertspace,findLZ,get_Lz_sphere
end

module TwoBodyGet
using WignerSymbols
function pp_matrix(s::Float64, m::Integer)
	dim = Int(2s+1)
	mat = zeros(dim,dim)
	for i in 1:dim
		for j in 1:dim
			#println((i,j))
			if abs(i+j-2-2s) <= (2s-m)
				mat[i,j] = clebschgordan(Float64,s,-i+s+1,s,-j+s+1,2s-m,-i-j+2s+2)
			end
		end
	end
    #@show mat
	return mat
end

function get_twobody(n_orb, v_list::Vector{T} where T<:Integer, c_list::Vector{Float64})
    twoBody = Vector{Tuple{Float64, Tuple{UInt64, UInt64,UInt64,UInt64}}}()
    s = (n_orb-1)/2
    vmat = [pp_matrix(s,v_list[i]) * sqrt(c_list[i]) for i in 1:length(v_list)]
    for i = 1:n_orb
        for j = 1:n_orb
            for l = 1:n_orb
                k = i+j-l
                if k>0 && k<n_orb+1
                    h = 0
                    for v in vmat
                        h += v[j,i]*v[l,k]
                    end
                    push!(twoBody,(h,(i,j,l,k)))
                end
            end
        end
    end
    return twoBody
end
export pp_matrix,get_twobody
end

module OneBodyGet
single_particle_state_sphere(θ::Number, φ::Number, S::Number, m::Number) = cos(θ/2)^(S+m) * sin(θ/2)^(S-m) * exp(m*im*φ) / sphere_coef(S,m)

single_particle_state_sphere(θ::Array{T} where T<: Number, φ::Array{T} where T<: Number, S::Number, m::Number) = cos.(θ./2).^(S.+m) .* sin.(θ/2).^(S.-m) .* exp.(m*im.*φ) / sphere_coef(S,m)

# Normalization coefficient on the sphere
sphere_coef(S,m) =   sqfactorial(S-m)/sqfactorial(S+m+1, 2S+1)


sqfactorial(N) = prod(map(sqrt, 1:N)) 
# square root of N!, avoid overflow for large N and more efficient than sqrt(factorial(big(N)))
# (overflow starts at N=21)

sqfactorial(n,N) = prod(map(sqrt, n:N))
function get_oneBody(n_orb::Integer,θ_list::Vector{Float64}, ϕ_list::Vector{Float64}, V_list::Vector{Float64})
    oneBody = Vector{Tuple{ComplexF64, Tuple{UInt64, UInt64}}}()
    S = (n_orb-1)/2
    for i = 1:length(θ_list)
        θ = θ_list[i]
        ϕ = ϕ_list[i]
        height = V_list[i]
        sfunction = map(m->single_particle_state_sphere(θ,ϕ,S,m), S:-1:-S)
        for j = 1:n_orb
            for k = 1:n_orb
                push!(oneBody,(conj(sfunction[j])*sfunction[k]*height,(j,k)))
            end
        end
    end
    return oneBody
end
using WignerD
function get_oneBody_widebump(n_orb::Int64,θ_list::Vector{Float64}, ϕ_list::Vector{Float64}, pinsize_list::Vector{Int64})
    oneBody = Vector{Tuple{ComplexF64, Tuple{UInt64, UInt64}}}()
    S = (n_orb-1)/2
    for i = 1:length(θ_list)
        θ = θ_list[i]
        ϕ = ϕ_list[i]
        W = WignerD.wignerD(S,0,θ,ϕ)
        pinsize = pinsize_list[i]
        for j = 1:n_orb
            for k = 1:n_orb
                for m0 = S:-1:S-pinsize+1
                    push!(oneBody,(W[Int(S-m0+1),j]*conj(W[Int(S-m0+1),k]),(j,k)))
                end
            end
        end
    end
    return oneBody
end

export get_oneBody,get_oneBody_widebump
end


module ConstructManybodyMatrix
using SparseArrays
function calcSign(I::UInt64, K::Integer)::Int64
    M = I-(I&(2K-1)) #occupied states to the left of (i)
    btwnCnt = count_ones(M) #counts all occupied to the left of (i)
    sign = (-1)^btwnCnt
    return sign
end

function calcSignNegative(I::UInt64, K::Integer)::Bool # Returns true if sign is negative, false if positive
    M = I-(I&(2K-1)) #occupied states to the left of (i)
    btwnCnt = count_ones(M) #counts all occupied to the left of (i)
    sign = (-1)^btwnCnt
    return sign < 0
end

#annihilation operator
#this function returns (result, sign, annhihilate)
#the second variable is false if sign is positive and true if sign is negative
#the third variable is false if the state is not annihilated and true if the state is annihilated
function C(i::Integer, I::UInt64)::Tuple{UInt64, Bool,Bool}
    #if sign == 0 #if sign is 0, the many body state is anihilated so can just proceed
    #    return 0, 0
    #end
    K = 2^(i-1)  #integer representation of single particle state
    if (I & K) != K # if (i) was already unoccupied, then anihilate many body state
        return I,true,true
    else
        L = I - K  #anihilate (i) in I if (i) is occupied in I
        sign = calcSignNegative(I, K)
        return L, sign,false
    end
end

#creation operator
function CDag(i::Integer, I::UInt64)::Tuple{UInt64, Bool,Bool}
    #if sign == 0 #if sign is 0, the many body state is anihilated so can just proceed
    #    return 0, 0
    #end
    K = 2^(i-1)  #integer representation of single particle state
    if I & K == K # if state is already occupied, anihilate many body state
        return I,true,true
    else
        L = I + K #create (i) in I
        sign = calcSignNegative(I, K)
        return L, sign,false
    end
end

#two-body part of Hamiltonian
function calcV(twoBody::Vector{Tuple{Float64, Tuple{UInt64, UInt64,UInt64, UInt64}}}, States::Vector{UInt64},subzone::Integer;num_zones=10,name="matrix")
    if !isdir("Matrix") mkdir("Matrix") end
    if !isdir("Matrix/two-body") mkdir("Matrix/two-body") end
    if !isdir("Matrix/two-body/$name") mkdir("Matrix/two-body/$name") end

    Rows = [Vector{UInt64}() for _ in 1:Threads.nthreads()]
    Cols = [Vector{UInt64}() for _ in 1:Threads.nthreads()]
    Vals = [Vector{ComplexF64}() for _ in 1:Threads.nthreads()]

    base_size,remainder = divrem(length(twoBody),num_zones)
    if subzone <= remainder
        begin_index = (base_size+1)*(subzone-1)+1
        end_index = (base_size+1)*subzone
    else
        begin_index = (base_size+1)*remainder + (subzone-remainder-1)*base_size+1
        end_index = (base_size+1)*remainder + (subzone-remainder)*base_size
    end
    nproc = end_index - begin_index + 1
    Threads.@threads for i in begin_index:end_index
        #print("\r$(i-begin_index+1)\t/\t$nproc         ")
        (v, (i, j, l, k)) = twoBody[i]
        for IInd in eachindex(States)
            I = States[IInd]
            totalsign = 0

            J, sign,annihilated = C(k, I)
            if annihilated continue end
            totalsign += sign
    
            J, sign, annihilated = C(l, J)
            if annihilated continue end
            totalsign += sign
            
            J, sign, annihilated = CDag(j, J)
            if annihilated continue end
            totalsign += sign

            J, sign,annihilated = CDag(i, J)
            if annihilated continue end
            totalsign += sign

            JInd = searchsortedfirst(States,J)
            if JInd == (length(States)+1) || States[JInd]!=J
                continue
            end

            #print("\r   $i\t$j\t$IInd\t$JInd\t\t")
            push!(Rows[Threads.threadid()-1], JInd)
            push!(Cols[Threads.threadid()-1], IInd)
            push!(Vals[Threads.threadid()-1], ((-1)^totalsign)*v) 
        end
    end
    rows=reduce(vcat,Rows)
    cols=reduce(vcat,Cols)
    vals=reduce(vcat,Vals)
    #FREE UP RAM
    Rows=nothing
    Cols=nothing
    Vals=nothing
    GC.gc() #garbage collect
    open("Matrix/two-body/$name/rows.txt","a+") do f
        write(f,join(string.(rows),"\n"))
        if subzone != num_zones
            write(f,"\n")
        end
    end
    open("Matrix/two-body/$name/cols.txt","a+") do f
        write(f,join(string.(cols),"\n"))
        if subzone != num_zones
            write(f,"\n")
        end
    end
    open("Matrix/two-body/$name/vals.txt","a+") do f
        write(f,join(string.(vals),"\n"))
        if subzone != num_zones
            write(f,"\n")
        end
    end
    
    rows=nothing
    cols=nothing
    vals=nothing
    GC.gc() 
    return
end
    
function calcT(oneBody::Vector{Tuple{ComplexF64, Tuple{UInt64, UInt64}}}, States::Vector{UInt64},subzone::Integer;num_zones=10,name="matrix")
    if !isdir("Matrix") mkdir("Matrix") end
    if !isdir("Matrix/one-body") mkdir("Matrix/one-body") end
    if !isdir("Matrix/one-body/$name") mkdir("Matrix/one-body/$name") end

    Rows = [Vector{UInt64}() for _ in 1:Threads.nthreads()]
    Cols = [Vector{UInt64}() for _ in 1:Threads.nthreads()]
    Vals = [Vector{ComplexF64}() for _ in 1:Threads.nthreads()]
    base_size,remainder = divrem(length(oneBody),num_zones)
    if subzone <= remainder
        begin_index = (base_size+1)*(subzone-1)+1
        end_index = (base_size+1)*subzone
    else
        begin_index = (base_size+1)*remainder + (subzone-remainder-1)*base_size+1
        end_index = (base_size+1)*remainder + (subzone-remainder)*base_size
    end
    nproc = end_index - begin_index + 1
    Threads.@threads for i in begin_index:end_index
        #print("\r$(i-begin_index+1)\t/\t$nproc         ")
        (t, (i, j)) = oneBody[i]
        for IInd in eachindex(States)
            I = States[IInd]
            totalsign = 0
            
            J, sign, annihilated = C(j, I)
            if annihilated continue end
            totalsign += sign

            J, sign, annihilated = CDag(i, J)
            if annihilated continue end
            totalsign += sign

            JInd = searchsortedfirst(States,J) #REQUIRES THAT States BE SORTED. see http://www.jlhub.com/julia/manual/en/function/searchsortedfirst
            if JInd == (length(States)+1) || States[JInd]!=J
                continue
            end
            push!(Rows[Threads.threadid()-1], JInd)
            push!(Cols[Threads.threadid()-1], IInd)
            push!(Vals[Threads.threadid()-1], ((-1)^totalsign)*t) 
        end
    end
    rows=reduce(vcat,Rows)
    cols=reduce(vcat,Cols)
    vals=reduce(vcat,Vals)
    #FREE UP RAM
    Rows=nothing
    Cols=nothing
    Vals=nothing
    GC.gc() #garbage collect
    
    open("Matrix/one-body/$name/rows.txt","a+") do f
        write(f,join(string.(rows),"\n"))
        if subzone != num_zones
            write(f,"\n")
        end
    end
    open("Matrix/one-body/$name/cols.txt","a+") do f
        write(f,join(string.(cols),"\n"))
        if subzone != num_zones
            write(f,"\n")
        end
    end
    open("Matrix/one-body/$name/vals.txt","a+") do f
        write(f,join(string.(vals),"\n"))
        if subzone != num_zones
            write(f,"\n")
        end
    end

    rows=nothing
    cols=nothing
    vals=nothing
    GC.gc() 

    return
end


export calcV,calcT,CDag,C
end


