#Requires Markov_chain.jl!

@doc raw"""`δ` calculates the square of the stochastic estimate for the error on the stochastic estimate 
    for the expectation value ⟨X⟩ of a given observable X(s) with thermalisation time `N_th`. 
    The stochastic estimate of ⟨X⟩ is the mean X̅ of X (doc of function `covariance_method`)
    calculated with N configurations in the Markov chain. For N → ∞ the stochastic estimate X̅ = ⟨X⟩. 
    If N is finite one can estimate the error of X̅ by the square root of `δ`(W):
    ```math
    \delta(W) = \frac{1}{(N-N_{th})^2}\sum_{N_{th} < n,m \le N\, \cap\, |m-n|\le W}[X(s_m)-\bar{X}_N][X(s_n)-\bar{X}_N]
    ```
    whereby W is the summation window.
    ##### Input parameters of `δ`:
    - `X::Array` : Observable X
    - `N_th::Int64` : Thermalisation time of observable X
    - `W::Int64`  : Some integrer used as summation window
    ##### Output values: 
    - `δ::Number` : Square of the stochastic estimate for the error as discribed above""" ->
function δ(X::Vector, N_th::Int64, W::Int64)
    @assert(N_th > 0, "N_th = " * string(N_th) * "must be positive")
    @assert(N_th < length(X), "N_th = " * string(N_th) * "must be smaller than N_config")
    @assert(W > 0, "W = " * string(W) * "must be positive")
    
    N_config = length(X)
    X_bar = mean(X[N_th+1:end])
    δX = X[N_th+1:end] .- X_bar
    out = 0
    
    for n in eachindex(δX)
        m_1, m_2 = max(1, n-W), min(N_config-N_th, n+W)
        out += δX[n] * sum(δX[m_1:m_2])
    end
    
    return 1/(N_config-N_th)^2 * out
end

@doc raw"""`covariance_method` calculates the stochastic estimate for the expectation value ⟨X⟩ of a given observable X(s) with thermalisation time `N_th`. 
    It also returns the stochastic estimate of its error, calculated with function `δ`. 
    The stochastic estimate of ⟨X⟩ is the mean 
    ```math
    \bar{X}_N = \frac{1}{N-N_{th}}\sum_{n=N_{th}+1}^NX(s_n)
    ```
    whereby N > `N_th` is the number of configurations of the `Markov_chain`. For N → ∞ the stochastic estimate X̅ = ⟨X⟩. 
    If N is finite one can estimate the error of X̅ by the square root of `δ`(W), whereby W is the summation window. 
    This restriction can be applied, because the summands (expression of the covariance matrix) become very close to 
    zero if |n-m| is large. 
    This function `covariance_method` also chooses W=W̅ such that the error made by neglecting the terms with |m-n|>W̅ is small. For this instance it requires  
    ```math
    |\frac{\delta(\bar{W}+1)-\delta(\bar{W})}{\delta(\bar{W})}| < \epsilon
    ```  
    whereby ϵ is some tolerance value, here 5% if not chosen otherwise.
    ##### Input parameters of `covariance_method`:
    - `X::Array` : Observable X
    - `N_th::Int64` : Thermalisation time of observable X
    - `W::Int64`  : Some integrer used as summation window
    ##### Optional parameters:
    - `tol::Float64=0.05` : Tolerance value as discribed above
    - `verbose::Bool=false` : If true chosen value of W = W̅ will be printed
    ##### Output values: 
    - `X̅::Number`  : Stochastic estimate of ⟨X⟩
    - `√δ::Number` : Stochastic estimate for the error of X̅""" ->
function covariance_method(X::Array, N_th::Int64, W_max::Int64; tol::Float64=0.05, verbose::Bool=false)
    W = 1
    δ_1 = δ(X, N_th, W)
    δ_2 = δ(X, N_th, W+1)
    
    while abs((δ_2-δ_1)/δ_1) >= tol && W < W_max
        W += 1
        δ_1, δ_2 = δ(X, N_th, W), δ(X, N_th, W+1)
    end
    
    if abs((δ_2-δ_1)/δ_1) >= tol error("Desired tolerance not achieved") end
    
    if verbose println("W_bar = " * string(W)) end
    
    return mean(X[N_th+1:end]), sqrt(δ_1)
end