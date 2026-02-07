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

function check_valid_result(dir_name) # A directory must contains "out.log", otherwise its result is invalid
    return ("out.log" in readdir(dir_name))
end

function read_rows_and_cols(one_body_file="",two_body_file="")
    rows = UInt64[]
    cols = UInt64[]
    ut   = BitVector([])
    lens = [0,0]
    try # two-body is read first
        open("Matrix/two-body/$(two_body_file)/rows-cols.txt","r") do f
            for line in eachline(f)
                row, col = parse.(Int,split(line,","))
                push!(rows,row)
                push!(cols,col)
                lens[2] += 1
            end
        end
        open("Matrix/one-body/$(one_body_file)/rows-cols.txt","r") do f
            for line in eachline(f)
                row, col = parse.(Int,split(line,","))
                push!(rows,row)
                push!(cols,col)
                #if row != col
                #    push!(rows, col)
                #    push!(cols,row)
                #end
                lens[1] += 1
            end
        end
        return rows, cols,lens
    catch SystemError
        println("One or more specified files not found in:")
        println(" • Matrix/one-body/$one_body_file")
        println(" • Matrix/two-body/$two_body_file")
        return rows, cols,lens
    end
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
                for subzone = 1:10
                    print("\rProgress $subzone/10\t")
                    calcV(twoBody,basis,subzone;num_zones=10,name=latest_file_2bdy)
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
                for subzone = 1:100
                    print("\rProgress $subzone/100\t")
                    calcT(oneBody,basis,subzone;num_zones=100,name=latest_file_1bdy)
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
            rows, cols, lens = read_rows_and_cols(latest_file_1bdy,latest_file_2bdy)
            #println("Number of non-zero elements: $(lens[1]), $(lens[2])")
            if sum(lens)>0
                H_matrix   = spzeros(ComplexF64,rows,cols,d,d) # requires Julia 1.10 or newer
                # Read one-body matrix values
                shiftindex = lens[2]
                #println("Shift index = $shiftindex")
                open("Matrix/one-body/$(latest_file_1bdy)/vals.txt") do f
                    for (i,line) in enumerate(eachline(f))
                        linevalue = parse(ComplexF64,line)
                        H_matrix[rows[i+shiftindex],cols[i+shiftindex]] += λ*linevalue
                        #if rows[i+shiftindex] != cols[i+shiftindex]
                        #    shiftindex += 1
                            #println("Updated shift index = $shiftindex")
                        #    H_matrix[rows[i+shiftindex],cols[i+shiftindex]] += λ*conj(linevalue)
                        #end
                        #H_matrix[cols[i+lens[2]],rows[i+lens[2]]] += λ*conj(linevalue)
                    end
                end
                # Read two-body matrix values
                open("Matrix/two-body/$(latest_file_2bdy)/vals.txt") do f
                    for (i,line) in enumerate(eachline(f))
                        H_matrix[rows[i],cols[i]] += parse(ComplexF64,line)
                    end
                end
            end
        end # end of @time block 

        rows = nothing
        cols = nothing
        GC.gc()
        
        if d<20
            display(real.(Matrix(H_matrix)))
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
