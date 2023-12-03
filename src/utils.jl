# https://discourse.julialang.org/t/smoothing-noisy-data-using-moving-mean/65329/6
function moving_average(A::AbstractArray, m::Int)
    out = similar(A)
    R = CartesianIndices(A)
    Ifirst, Ilast = first(R), last(R)
    I1 = m÷2 * oneunit(Ifirst)
    for I in R
        n, s = 0, zero(eltype(out))
        for J in max(Ifirst, I-I1):min(Ilast, I+I1)
            s += A[J]
            n += 1
        end
        out[I] = s/n
    end
    return out
end

"""
    hampel(x, k, [,thr])

    Simple Hampel Filter for Julia lang.

    # Arguments
    - `x::n-elements Vector{Float64}`: time series data.
    - `k::Int`: the step size of a filter.
    - `thr::Float64`: Float threshold value (default: 3)

    # Examples

    ```jldoctest
        julia> x = [1, 2, 3, 10, 3, 2];
        julia> hampel(x,3)
        6-element Vector{Float64}:
        1.0
        2.0
        3.0
        2.5
        3.0
        2.0
    ```
"""

function hampel(x,k,thr=3)
    
    # Exceptions
    if length(x)==0
        throw(DomainError(x, "hampel() is not defined for zero-length vector"))
    elseif ~isinteger(k) || k<1
        throw(DomainError(k, "k should be positive integer."))
    elseif thr<=0
        throw(DomainError(thr, "thr should be positive."))
    end

    arrSize = length(x)
    output_x = Float64[]
    for i in 1:arrSize
        if i<2k
            kernel = collect(1:i)
        elseif i+2k>arrSize
            kernel = collect(i:arrSize)
        else
            kernel = collect(i:i+2k)
        end
        med = median(x[kernel])
        std = 1.4826*median(abs.(x[kernel].-med))
        if abs(x[i]-med)>thr*std
            push!(output_x,med)
        else
            push!(output_x,x[i])
        end
    end
    output_x
end