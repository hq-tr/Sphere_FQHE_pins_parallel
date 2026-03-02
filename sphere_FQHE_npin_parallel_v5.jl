include("SphereED_v5.jl")
include("FQH_states_IO.jl")
using .FQH_states
using .HilbertSpaceGenerator
using .OneBodyGet
using .TwoBodyGet
using .ConstructManybodyMatrix
using LinearAlgebra
using SparseArrays
using ArgMacros
using BenchmarkTools
using Printf
using WignerSymbols
using JLD2
using Arpack
using Dates
using Combinatorics


function check_valid_result(dir_name) # A directory must contains "out.log", otherwise its result is invalid
    return ("out.log" in readdir(dir_name))
end

dec2dex(n::Integer) = [i-1 for (i,d) in enumerate(digits(n,base=2)) if d==1]

function get_rows_and_cols(basis::Vector{T} where T<:Integer,n_orb::Integer;onebody=true,twobody=true)
    if !onebody && !twobody return UInt64[],UInt64[] end

    d = length(basis)

    #for b in basis
    #   println(bitstring(b))
    #end

    # Initialize the rows and cols vectors with diagonal terms
    rows = UInt64.(1:d)
    cols = UInt64.(1:d)

    # For every element in basis, find all other elements with which it gives non-zero matrix elements
    for (index1,basis1) in enumerate(basis)
        #print("\r   Checking basis $(index1) of $d\t\t\t")
        basis1_dex = dec2dex(basis1)

        # One-body c†_m2 c_m1
        if onebody
            for m1 in basis1_dex
                for m2 in (m1+1):(n_orb-1) #only check m2 > m1
                    
                    if m2 in basis1_dex continue end # Make sure m2 is not an electron

                    basis2 = basis1 - 2^m1 + 2^m2
                    index2 = searchsortedfirst(basis,basis2)

                    push!(rows,index1)
                    push!(cols,index2)
                    push!(rows,index2)
                    push!(cols,index1)

                end
            end
        end

        # Two-body c†_m3 c†_m4 c_m1 c_m2 (m3 < m1 < m2 < m4)
        if twobody
            for (m1,m2) in combinations(basis1_dex,2)
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
    end
    return rows,cols
end

function main()
    # ================================ READ USER INPUT ================================
    @inlinearguments begin
        @argumentdefault Int 5 k "-n" "--nev"
        @argumentdefault String "" intname "-i" "--interaction-file"
        @argumentdefault Int 1 npin "--npins"
        @argumentrequired Int n_el "-e" "--n_el"
        @argumentrequired Int n_orb "-o" "--n_orb"
        @argumentdefault Int 1 pin_size "-k" "--pin-size"
        @argumentdefault Float64 1.0 λ "--lambda"
        @argumentflag forcerebuild "--force-rebuild-matrix"
        @argumentflag nooutput "--no-output"
        @argumentflag noeigenstate "--no-eigenstate"
        @argumentflag nomatrix "--no-saving-matrix"
        @argumentdefault Int 1 numzones1bd "-z" "--num-zones-1bdy"
        @argumentdefault Int 1 numzones2bd "-Z" "--num-zones-2bdy"
    end
    #LzConserve = npin ≤ 2
    LzConserve = false # This version doesn't work for rotationally symmetric configurations yet
    #λ = lamda_list[nlamda]
    println("============================================================")
    println("      FULL-ED OF TWO-BODY INTERACTION ON THE SPHERE")
    println("               with potential pins")
    println()
    println("This version loads the existing two-body matrix from file (if it exists)")
    println("and build the one-body matrix on top of it.")
    println("============================================================")


    # Get the current time
    timenow = now()
    println("Local time and date: $timenow")
    # Reading basis input
    if n_el != nothing && n_orb != nothing
        if LzConserve
            println("Generating a basis with $(n_el) electrons and $(n_orb) orbitals with different Lz sectors.")
            basis = Vector{Vector{Int64}}()
            Lz_max = sum((n_orb-1)/2:-1:((n_orb-1)/2-n_el+1))
            for Lz in Lz_max:-1:-Lz_max
                push!(basis,fullhilbertspace(n_el,n_orb,Lz))
            end
        else
            println("Generating a basis with $(n_el) electrons and $(n_orb) orbitals (all Lz sectors).")
            basis = fullhilbertspace(n_el,n_orb)
            println("The dimension is $(length(basis))")
        end
        outname = @sprintf "%ie_%io" n_el n_orb
    else
        println()
        println("WARNING: No input or incomplete input was specified. The program will now terminating.")
        println("Run the program with '-h' or '--help' tag to view possible arguments.")
        println()
        return
    end

    # Reading two-body interaction input
    v_list = Int32[]
    c_list = Float64[]

    if intname != "none"
        if length(intname) == 0
            println("Input m for Vₘ and the corresponding coefficient. ")
            println("Each pp term takes one line, with two numbers separated by a space.")
            println("Put a 0 to end")
            reading = true
            while reading
                data = readline()
                if data == "0"
                    reading = false
                else
                    try
                        pp = split(data)
                        push!(v_list,parse(Int32, pp[1]))
                        push!(c_list,parse(Float64,pp[2]))
                    catch
                        println("Invalid input. Try again or input 0 to end.")
                    end
                end
            end
        else
            println("Reading interaction from $(intname).")
            if isfile(intname)
                open(intname) do f
                    for line in map(s->split(s),readlines(f))
                        append!(v_list,parse(Int32,line[1]))
                        append!(c_list,parse(Float64,line[2]))
                    end
                end
            else
                print("Interaction file '$(intname)' not found. Terminating.")
                return false
            end
        end
    end

    #set pins
    if npin == 1
        θ_list = [0.0]
        ϕ_list = [0.0]
    elseif npin == 2
        θ_list = [0.0,π]
        ϕ_list = [0.0,0.0]
    elseif npin == 3
        θ_list = [0.0,2/3*π,4/3*pi]
        ϕ_list = [0.0,0.0,0.0]
    elseif npin == 4
        θ₁ = π/2 + atan(1/√(8))
        ϕ₁ = 0.
        ϕ₂ = 2π/3
        ϕ₃ = 4π/3
        θ_list = [0,θ₁,θ₁,θ₁]
        ϕ_list = [0,ϕ₁,ϕ₂,ϕ₃]
    elseif npin == 6
        θ_list = [0,π/2,π/2,π/2,π/2,π]
        ϕ_list = [0,0,π/2,π,3π/2,0]
    end
    pinsize_list = ones(Int,npin)*pin_size

    # ======================== CONSTRUCT AND DIAGONALIZE HAMILTONIAN ======================
    println("--------")
    println("Constructing the Hamiltonian")
    outname = @sprintf "%ie_%io" n_el n_orb
    dirname = "wide_pinnumber_$(length(pinsize_list))_pinsize_$(pin_size)_out"
    if !isdir(dirname) mkdir(dirname) end

    oneBody = get_oneBody_widebump(n_orb,θ_list, ϕ_list, pinsize_list)
    twoBody = get_twobody(n_orb, v_list::Vector{Int32}, c_list::Vector{Float64})
    if LzConserve
	   for nLz = 1:length(basis)
            d = length(basis[nLz])
            if !SaveMemory
                H_2bdy = calcV(twoBody,basis[nLz])
            else
                H_2bdy = spzeros(ComplexF64, d,d);
                for subzone = 1:10
                    print("\rProgress $subzone/10\t")
                    H_2bdy = calcV(H_2bdy,twoBody,basis[nLz],subzone);
                end
            end
            H_1bdy = calcT(oneBody,basis[nLz])
            for λ in lambda_list
                H_matrix = H_2bdy+ λ*H_1bdy
                if d<10
                    H_matrix = Matrix(H_matrix)
                    ϵ, ϕ = eigen(H_matrix)
                else
                    ϵ, ϕ = eigs(H_matrix, nev=k,which=:SR)
                end
                open("$(dirname)/eigen_$(outname)_λ_$λ.txt","a+") do f
                    for i in 1:min(k,d)
                        gs_coef =  ϕ[:,i] 
                        write(f,"$(ϵ[i])\t")
                        LZ = get_Lz_sphere(basis[nLz],gs_coef,n_orb)
                        write(f,"$(LZ)\n")
                    end
                end
            end
        end
    else
        d   = length(basis) # Dimension
        println(" • Two-body")
        build_matrix = true
        latest_file_2bdy = "$(n_el)e_$(n_orb)o_$(intname)_$timenow"
        if !forcerebuild
            if isdir("Matrix")
                if isdir("Matrix/two-body")
                    matrix_directories = filter(x -> occursin("$(n_el)e_$(n_orb)o_$intname",x),readdir("Matrix/two-body"))
                    while length(matrix_directories) > 0
                        if check_valid_result("Matrix/two-body/$(matrix_directories[end])")
                            latest_file_2bdy = matrix_directories[end]
                            build_matrix = false
                            break
                        else
                            pop!(matrix_directories)
                        end
                    end
                else
                    build_matrix = true
                end
            else 
                build_matrix = true
            end
        end
        
        if build_matrix
            @time begin
                for subzone = 1:numzones2bd
                    print("\rProgress $subzone/$numzones2bd\t")
                    calcV(twoBody,basis,subzone;num_zones=numzones2bd,name=latest_file_2bdy)
                end
            end 
            open("Matrix/two-body/$latest_file_2bdy/out.log","w+") do f
                write(f,"Completed on $(now()).")
            end
        else
            println("Matrix data found at")
            println("  Matrix/two-body/$latest_file_2bdy")
        end

        println(" • One-body")
        build_matrix = true
        latest_file_1bdy = "$(n_el)e_$(n_orb)o_$(npin)_pins_k_$(pin_size)_$timenow"
        if !forcerebuild
            if isdir("Matrix")
                if isdir("Matrix/one-body")
                    matrix_directories = filter(x -> occursin("$(n_el)e_$(n_orb)o_$(npin)_pins_k_$(pin_size)",x),readdir("Matrix/one-body"))
                    while length(matrix_directories) > 0
                        if check_valid_result("Matrix/one-body/$(matrix_directories[end])")
                            latest_file_1bdy = matrix_directories[end]
                            build_matrix = false
                            break
                        else
                            pop!(matrix_directories)
                        end
                    end
                else
                    build_matrix = true
                end
            else 
                build_matrix = true
            end
        end


        if build_matrix
            @time begin
                for subzone = 1:numzones1bd
                    print("\rProgress $subzone/$(numzones1bd)\t")
                    calcT(oneBody,basis,subzone;num_zones=numzones1bd,name=latest_file_1bdy)
                end
            end      
            open("Matrix/one-body/$latest_file_1bdy/out.log","w+") do f
                write(f,"Completed on $(now()).")
            end
        else
            println("Matrix data found at")
            println("  Matrix/one-body/$latest_file_1bdy")
        end

        println("--------------------")
        println("Reading matrix from files")
        @time begin
            # Construct the matrix from saved data
            rows, cols = get_rows_and_cols(basis,n_orb)
            nnz = length(rows)
            

            H_matrix   = spzeros(ComplexF64,rows,cols,d,d) # requires Julia 1.10 or newer
            
            rows = nothing
            cols = nothing
            GC.gc() # free up memory

            println("A sparse matrix was pre-allocated with $nnz non-zero elements")

            println("Begin reading matrix elements")

            # Read one-body matrix values
            f = open("Matrix/one-body/$(latest_file_1bdy)/rows-cols.txt")
            g = open("Matrix/one-body/$(latest_file_1bdy)/vals.txt")
            
            for (i,(linef,lineg)) in enumerate(zip(eachline(f),eachline(g)))
                print("\r  Reading line $i of one-body matrix\t\t")
                if length(linef) > 0 && length(lineg) > 0
                    row, col = parse.(UInt64,split(linef,","))
                    linevalue = parse(ComplexF64,lineg)

                    H_matrix[row,col] += λ*linevalue
                    H_matrix[col,row] += λ*conj(linevalue)
                end
            end

            close(f)
            close(g)
            println("Done")

            # Read two-body matrix values
            f = open("Matrix/two-body/$(latest_file_2bdy)/rows-cols.txt")
            g = open("Matrix/two-body/$(latest_file_2bdy)/vals.txt")
            
            for (i,(linef,lineg)) in enumerate(zip(eachline(f),eachline(g)))
                print("\r  Reading line $i of two-body matrix\t\t")
                if length(linef) > 0 && length(lineg) > 0
                    row, col = parse.(UInt64,split(linef,","))
                    H_matrix[row,col] += parse(ComplexF64,lineg)
                end
            end

            close(f)
            close(g)
            println("Done")
        end # end of @time block 
        
        if d<20
            display(Matrix(H_matrix))
        elseif d<500
            display(H_matrix)
        else
            println(summary(H_matrix))
        end

        println("--------------------")
        println("Diagonalizing with ARPACK")

        ϵ, ϕ = eigs(H_matrix, nev=k,which=:SR)
        println("Eigenvalues = ")
        for ee in abs.(ϵ)
            println("  $ee")
        end

        if !nooutput
            open("$(dirname)/eigen_$(outname)_λ_$λ.txt","w+") do f
               #write(f,"$λ\n")
               for i in 1:k
                   gs_coef =  ϕ[:,i] 
                   write(f,"$(abs(ϵ[i]))\n")
                   #LZ = get_Lz_sphere(basis,gs_coef,n_orb)
                   #write(f,"$(LZ)\n")
               end
            end
            if !noeigenstate
                for i in 1:k
                    gs = FQH_state(basis,ϕ[:,i],n_orb)
                    printwf(gs;fname="$(dirname)/gs_$(outname)_λ_$(λ)_$(i-1)")
                end
            end
        end

        # Delete the matrix files if required
        if nomatrix
            rm("Matrix/two-body/$(latest_file_2bdy)",recursive=true)
            rm("Matrix/one-body/$(latest_file_1bdy)",recursive=true)
        end
    end
end

@time main()
