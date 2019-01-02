"""
    sri(drift::Function, diffusion::Function,
        𝓘::AbstractRange,
        m::Integer, X₀::AbstractVector{<:Real},
        dW::Union{AbstractArray,Nothing}=nothing)

TODO: write docs

"""
function sri(drift::Function, diffusion::Function,
             𝓘::AbstractRange,
             m::Integer, X₀::AbstractVector{<:Real},
             dW::Union{AbstractArray,Nothing}=nothing)

    h = step(𝓘)
    d = length(X₀)
    common_type = promote_type(eltype(X₀),Float64) # get at least Float64
    if dW == nothing
        _dW = similar(X₀,common_type,m)
    elseif size(dW) ≠ (length(𝓘),m)
        error("Number of time points and Brownian increments do not match!")
    end

    # Preallocate
    Y = similar(X₀,common_type,length(𝓘),d)
    Y[1,:] = X₀
    I = similar(X₀,common_type,m,m)
    aH₁ = similar(X₀,common_type,d)
    bH₁ = similar(X₀,common_type,d,m)
    bH₂ = similar(X₀,common_type,d,m)
    bH₃ = similar(X₀,common_type,d,m)


    for n = 1:length(𝓘)-1
        if dW == nothing
            randn!(_dW)
            lmul!(√h, _dW)
        else
            _dW = view(dW,n+1,:)
        end
        I = simdoubleintegrals(_dW, h) / √h


        aH₁ = drift(𝓘[n], Y[n,:])
        bH₁ = diffusion(𝓘[n], Y[n,:])
        H₂(k) = Y[n,:] + bH₁ * I[:,k]
        H₃(k) = Y[n,:] - bH₁ * I[:,k]
        Y[n+1,:] .= Y[n,:] .+ aH₁.*h .+ bH₁*_dW
        for k = 1:m
            bH₂ = diffusion(𝓘[n], H₂(k))
            bH₃ = diffusion(𝓘[n], H₃(k))
            Y[n+1,:] .+= 0.5.*√h.* ( bH₂[:,k] .- bH₃[:,k] )
        end

        if any(isinf.(Y[n+1,:]))
            println("Got Inf! Stopping now...")
            Y[n+2:end,:] .= Y[n+1,:]'
            # @show n length(𝓘) Y[n-5:n,:] aH₁ h bH₁ _dW Y[n+1,:]
            return Y
        end
    end
    Y
end
"""
    sri for one-dimensional processes
    this needs out-of-place drift & diffusion
"""
function sri(drift::Function, diffusion::Function,
             𝓘::AbstractRange,
             m::Integer, X₀::Real,
             dW::Union{AbstractArray,Nothing}=nothing)

    h = step(𝓘)
    common_type = promote_type(typeof(X₀),Float64) # get at least Float64
    if dW == nothing
        _dW = Vector{common_type}(undef,m)
    elseif size(dW) ≠ (length(𝓘),m)
        error("Number of time points and Brownian increments do not match!")
    end

    # Preallocate
    Y = Vector{common_type}(undef,length(𝓘))
    Y[1] = X₀
    I = Matrix{common_type}(undef,m,m)
    bH₁ = Vector{common_type}(undef,m)
    bH₂ = Vector{common_type}(undef,m)
    bH₃ = Vector{common_type}(undef,m)

    for n = 1:length(𝓘)-1
        if dW == nothing
            randn!(_dW)
            lmul!(√h, _dW)
        else
            _dW = view(dW,n+1,:)
        end
        I = simdoubleintegrals(_dW, h) / √h

        aH₁ = drift(𝓘[n], Y[n])
        bH₁ = diffusion(𝓘[n], Y[n])
        H₂(k) = Y[n] + bH₁' * I[:,k]
        H₃(k) = Y[n] - bH₁' * I[:,k]
        Y[n+1] = Y[n] + aH₁.*h + bH₁'*_dW
        for k = 1:m
            bH₂ = diffusion(𝓘[n], H₂(k))
            bH₃ = diffusion(𝓘[n], H₃(k))
            Y[n+1] += 0.5*√h* ( bH₂[k] - bH₃[k] )
        end

        if isinf(Y[n+1])
            println("Got Inf! Stopping now...")
            Y[n+2:end] .= Y[n+1]
            # @show n length(𝓘) Y[n-5:n,:] aH₁ h bH₁ _dW Y[n+1,:]
            return Y
        end
    end
    Y
end

function sri_inplace(drift::Function, diffusion::Function,
             𝓘::AbstractRange,
             m::Integer, X₀::AbstractVector{<:Real})

    h = step(𝓘)
    d = length(X₀)

    # Preallocate
    Y = zeros(length(𝓘),d)
    Y[1,:] = X₀
    I = zeros(m,m)
    aH₁ = zeros(d)
    bH₁ = zeros(d,m)
    bH₂ = zeros(d,m)
    bH₃ = zeros(d,m)

    for n = 1:length(𝓘)-1
        W = randn(m)
        if m > 1
            n_approx::Int64 = ceil(Int64, √( m*(m-1)*(m+4*(W'*W))/(h*24*π^2) ))
            I = simdoubleintegrals(W,n_approx)
        else
            I = simdoubleintegrals(W)
        end
        lmul!(√h,W)
        lmul!(h,I)

        drift(aH₁, 𝓘[n], Y[n,:])
        diffusion(bH₁, 𝓘[n], Y[n,:])
        H₂(k) = Y[n,:] + bH₁ * I[:,k]
        H₃(k) = Y[n,:] - bH₁ * I[:,k]
        Y[n+1,:] .= Y[n,:] .+ aH₁.*h .+ bH₁*W
        for k = 1:m
            diffusion(bH₂, 𝓘[n], H₂(k))
            diffusion(bH₃, 𝓘[n], H₃(k))
            Y[n+1,:] .+= 0.5.*√h.* ( bH₂[:,k] .- bH₃[:,k] )
        end
    end
    Y
end
