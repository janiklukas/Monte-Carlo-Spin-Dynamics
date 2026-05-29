#Requires LinearAlgebra,Random,Statistics,Dates,Printf!

@doc """The function `s_init` constructs the initial spin configuration on a `D`-dimensional lattice with `N` spins in each direction.
    ##### Input parameters of `s_init`:
    - `N::Int64`    : Number of spins in each direction
    - `D::Int64`    : Dimension of the lattice
    - `type::String` : Must be one of these
        - `"cold"` :  All spins equal to +1
        - `"hot"`  :  Randomly generated spin configuration (each spin has 50% probability to be +1 or -1)
    ##### Output values: 
    - `s::Array`  : D-dim Array with N Spins of +1 or -1 in each direction""" ->
function s_init(N::Int64, D::Int64, type::String)
    types = ["cold", "hot"]
    @assert(type in types, "type must be one of " * string(types))
    
    if type == "cold" return fill(1, Tuple(N for k in 1:D))
    else return rand([-1, 1], Tuple(N for k in 1:D)) end
end

@doc raw"""The function `Hamiltonian` calculates the energy of a spin configuration.
    ```math
    H(s) = -J\sum_n \sum_{k=1}^D s(\underline{n})s(\underline{n}+\underline{e}_k) - B\mu\sum_n s(\underline{n})
    ```
    whereby the first term corresponds to the nearest neighbor interaction and the second to the interaction with an external magnetic field b.
    ##### Input parameters of `Hamiltonian`:
    - `β::Float64` : β = J/kT dimensionless coupling constant
    - `b::Float64` : b = μB/kT dimensionless external magnetic field
    - `s::Array`  : D-dim Array with N Spins of +1 or -1 in each direction
    ##### Output values: 
    - `H::Float64` : Energy value as described above
    """ ->
function Hamiltonian(β::Float64, b::Float64, s::Array{Int64})
    s_plus = sum([circshift(s, 1*Matrix(I, D, D)[k,:]) for k in 1:D])
    return -β * sum(s .* s_plus) - b * sum(s)
end

@doc """The `Metropolis` algorithm constructs a new spin configuration. When applied multiple times it can be used in a `Markov_chain` (see function below).
    ##### Input parameters of `Metropolis`:
    - `β::Float64` : β = J/kT dimensionless coupling constant
    - `b::Float64` : b = μB/kT dimensionless external magnetic field
    - `s::Array{Int64}`  : D-dim Array with N Spins of +1 or -1 in each direction
    ##### Output values: 
    - `s::Array{Int64}`  : D-dim Array with new configuration of the N Spins in each direction""" ->
function Metropolis(β::Float64, b::Float64, s::Array{Int64})
    D = ndims(s)
    N = size(s)[1]
    
    #for each grid pt n, calculate the total spin of all pts of the form (n ± e_k) for some k:
    s_plus = sum([circshift(s, 1*Matrix(I, D, D)[k,:]) for k in 1:D])
    s_minus = sum([circshift(s, -1*Matrix(I, D, D)[k,:]) for k in 1:D])
    
    for i in CartesianIndices(s)
        ΔH = β*s[i] * (s_plus[i] + s_minus[i]) + b*s[i]
        
        if exp(-2*ΔH) > rand(Float64)
            s[i] *= -1
            for k in 1:D #update s_plus & s_minus                
                i_plus = ntuple(j -> j==k ? mod1(i[j] + 1, N) : i[j], D)
                i_minus = ntuple(j -> j==k ? mod1(i[j] - 1, N) : i[j], D)
                s_minus[CartesianIndex(i_plus)] += 2*s[i]
                s_plus[CartesianIndex(i_minus)] += 2*s[i]
            end
        end
    end
    
    return s
end

@doc raw"""The `Markov_chain`- Monte Carlo algorithm generates a stochastic sequence of magnetisation and energy configurations by iterating over the `Metropolis` function.
    The thermalized configuration is distributed with the probability
    ```math
    P(s) = \frac{1}{Z}e^{-\frac{H(s)}{k_B T}}
    ```
    This function needs `s_init`, `Hamiltonian` and `Metropolis` to work.
    ##### Input parameters of `Markov_chain`:
    - `β::Float64` : β = J/kT dimensionless coupling constant
    - `b::Float64` : b = μB/kT dimensionless external magnetic field
    - `N::Int64`    : Number of particles in each direction
    - `D::Int64`    : Dimension of the lattice
    - `init::string` : Type of initial spin configuration. Must be one of:
        - `"cold"` :  All spins equal to +1
        - `"hot"`  :  Randomly generated spin configuration (each spin has 50% probability to be +1 or -1)
    - `N_config::Int`: Number of configurations
    ##### Output values: 
    - `M::Vector{Int64}`  : Vector of length `N_config` containing the magnetisation for each configuration
    - `E::Vector{Float64}`  : Vector of length `N_config` containing the energy for each configuration""" ->
function Markov_chain(β::Float64, b::Float64, N::Int64, D::Int64, init::String, N_config::Int64)
    @assert(β > 0, "β = " *string(β) * " must be positive")
    @assert(b > 0, "b = " *string(b) * " must be positive")
    @assert(N_config > 0, "N_config = " * string(N_config) * " must be positive")

    s = s_init(N, D, init)
    M, E = [sum(s)], [Hamiltonian(β, b, s)]
    
    for m in 1:N_config
        s = Metropolis(β, b, s)
        M = vcat(M, sum(s))
        E = vcat(E, Hamiltonian(β, b, s))
    end
    
    return M, E
end

function plot_history(N_max::Int64, M::Vector, E::Vector; skip::Int64=10)
    x = collect(0:skip+1:N_max)
    p = plot(legend=:topright)

    scatter!(x, M[1:skip+1:N_max+1], color = 1, label = "M(s)")
    scatter!(x, E[1:skip+1:N_max+1], color = 2, label = "E(s)")
    xlabel!("Configuration number")

    return p
end