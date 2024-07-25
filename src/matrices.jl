module Matrix

# Transpose of a square matrix
function mT(m1::AbstractMatrix{T}) where T
    """
    Computes the transpose of a square matrix
    :param m1: rotation matrix
    :return mt: m1 transpose
    """
    n = size(m1, 1)
    m = size(m1, 2)
    mt = zeros(T, m, n)
    for ii in 1:n
        for jj in 1:m
            mt[jj, ii] = m1[ii, jj]
        end
    end
    return mt
end

# Matrix addition
function mxadd(m2::AbstractMatrix{T}, m1::AbstractMatrix{T}) where T
    return m2 .+ m1
end

# Matrix subtraction
function mxsub(m2::AbstractMatrix{T}, m1::AbstractMatrix{T}) where T
    return m2 .- m1
end

# Matrix multiplication; general case
function mxm(m2::AbstractMatrix{T}, m1::AbstractMatrix{T}) where T
    return m2 * m1
end

# Matrix multiplication for vector by a matrix
function mxv(m1::AbstractMatrix{T}, v1::AbstractVector{T}) where T
    return m1 * v1
end

# Scaling a matrix by a scalar
function mxs(scalar::T, m1::AbstractMatrix{T}) where T
    return scalar * m1
end

# Generates a skewed cross-product matrix from a vector
function skew(v1::AbstractVector{T}) where T
    @assert length(v1) == 3 "Input vector must be of length 3"
    v_tilde = zeros(T, 3, 3)
    v_tilde[1, 2] = -v1[3]
    v_tilde[1, 3] =  v1[2]
    v_tilde[2, 1] =  v1[3]
    v_tilde[2, 3] = -v1[1]
    v_tilde[3, 1] = -v1[2]
    v_tilde[3, 2] =  v1[1]
    return v_tilde
end

end # module
