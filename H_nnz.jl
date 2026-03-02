include("SphereED_v5.jl")
using .HilbertSpaceGenerator
using SparseArrays
using Combinatorics
using ArgMacros

dec2dex(n::Integer) = [i-1 for (i,d) in enumerate(digits(n,base=2)) if d==1]

function get_rows_and_cols(basis::Vector{T} where T<:Integer,n_orb::Integer)
    d = length(basis)

    for b in basis
       println(bitstring(b))
    end

    # Initialize the rows and cols vectors with diagonal terms
    rows = UInt64.(1:d)
    cols = UInt64.(1:d)

    # For every element in basis, find all other elements with which it gives non-zero matrix elements
    for (index1,basis1) in enumerate(basis)
        print("\r   Checking basis $(index1) of $d\t\t\t")
        basis1_dex = dec2dex(basis1)

        # One-body c†_m2 c_m1
        # for m1 in basis1_dex
        #     for m2 in (m1+1):(n_orb-1) #only check m2 > m1
                
        #         if m2 in basis1_dex continue end # Make sure m2 is not an electron

        #         basis2 = basis1 - 2^m1 + 2^m2
        #         index2 = searchsortedfirst(basis,basis2)

        #         push!(rows,index1)
        #         push!(cols,index2)
        #         push!(rows,index2)
        #         push!(cols,index1)

        #     end
        # end

        # Two-body c†_m3 c†_m4 c_m1 c_m2 (m3 < m1 < m2 < m4)
        for (m1,m2) in combinations(basis1_dex,2)
            #if m1+m2 ≥ n_orb continue end
            for m4 in (m2+1):min(n_orb,m2+m1)
                if m4 in basis1_dex continue end
                m3 = m1+m2-m4
                if m3 in basis1_dex continue end
                if m3 > n_orb continue end

                basis2 = basis1 - 2^m1 - 2^m2 + 2^m3 + 2^m4
                index2 = searchsortedfirst(basis,basis2)
                if index2 > d continue end
                if basis[index2] != basis2 continue end

                push!(rows,index1)
                push!(cols,index2)
                push!(rows,index2)
                push!(cols,index1)
            end
        end
    end
    return rows, cols
end

function main()
    println("\n-----------------------\n")
    println("This routine determines the locations of non-zero elements in 2-bdy + 1-bdy potential")
    println("\n-----------------------\n")

    @inlinearguments begin
    @argumentrequired Int n_el "-e" "--n_el"
    @argumentrequired Int n_orb "-o" "--n_orb"
    end

    println("Generating a basis with $(n_el) electrons and $(n_orb) orbitals (all Lz sectors).")
    @time basis = fullhilbertspace(n_el,n_orb)
    d = length(basis)
    println("The dimension is $d")

    @time rows,cols = get_rows_and_cols(basis,n_orb)
    println("Done")

    if d < 20
        for (row, col) in zip(rows,cols)
            println("$row\t$col")
        end
    end

    # Note that except for the diagonal terms, one-body and two-body potentials
    # should have different non-zero elements

    println("Matrix = ")
    H_mat = spzeros(ComplexF64,rows,cols,d,d)
    if d< 500
        display(H_mat)
    else
        println(summary(H_mat))
    end

end

@time main()