"""
    em(drift::Function, diffusion::Function,
        𝓘::AbstractRange,
        m::Integer, X₀::AbstractVector{<:Real},
        dW::Union{AbstractArray,Nothing}=nothing)

TODO: write docs

"""
function em(drift::Function, diffusion::Function,
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
            lmul!(√h,_dW)
        else
            _dW = view(dW,n+1,:)
        end

        aH₁ = drift(𝓘[n], Y[n,:])
        bH₁ = diffusion(𝓘[n], Y[n,:])
        Y[n+1,:] .= Y[n,:] .+ aH₁.*h .+ bH₁*_dW

        if any(isinf.(Y[n+1,:]))
            println("Got Inf! Stopping now...")
            Y[n+2:end,:] .= Y[n+1,:]'
            # @show n length(𝓘) Y[n-5:n,:] aH₁ h bH₁ _dW Y[n+1,:]
            return Y
        end
    end

    Y
end
