using LinearAlgebra

function adapt(f::Function, a::Real, b::Real)::Function
    return x -> f((b - a)*x/2.0 + (a + b)/2.0)*((b - a)/2.0)
end

function polyn(coeffs::AbstractVector, x::Real)::Real
    return sum(coeffs[i]*x^(i - 1) for i = 1:length(coeffs))
end

function polyndercoeffs(coeffs::AbstractVector)::AbstractVector
    return [(i-1)*coeffs[i] for i = 2:length(coeffs)]
end

function even(n::Integer)::Bool
    return n % 2 == 0
end

function legendrepolyncoeffs(n::Integer)::AbstractVector
    coeffs = zeros(n + 1)
    coeffs[1] = even(n) ? 1.0 : 0.0
    coeffs[2] = even(n) ? 0.0 : 1.0
 
    if n < 2
        return coeffs
    end

    i0 = even(n) ? 0 : 1
    
    for i = i0:n-2
        coeffs[i + 3] = -(n - i)*(n + i + 1)*coeffs[i + 1]/((i + 2)*(i + 1))
    end

    return coeffs ./ polyn(coeffs, 1.0)
end

function getcompanionmatrix(coeffs::AbstractVector)::AbstractMatrix
    n = length(coeffs)-1
    A = [zeros(1, n - 1) ; Matrix(I, n - 1, n - 1)]

    return [A -coeffs[1:n]./coeffs[n + 1]]
end

function getweight(n::Integer, x::Real)::Real
    coeffs = legendrepolyncoeffs(n)
    dercoeffs = polyndercoeffs(coeffs) 
    return 2.0/((1.0 - x^2)*(polyn(dercoeffs, x))^2)
end

function getpoints(n::Integer)::AbstractVector 
    coeffs = legendrepolyncoeffs(n)
    A = getcompanionmatrix(coeffs)
    return eigvals(A)
end

#NOTE: broadcast does a lot of unnecessary allocations, but it shouldn't cause
#too much performance issues
function getweights(n::Integer, points::AbstractVector)::AbstractVector
    return getweight.(n, points)
end

function quadsum(f::Function, points::AbstractVector, weights::AbstractVector)::Real
    return sum(f.(points).*weights)
end

function lgquad(f::Function, n::Integer, a::Real, b::Real)::Real
    points = getpoints(n)
    if any(typeof.(points) .<: Complex)
        error("points are complex")
    end

    weights = getweights(n, points)

    if a == -1.0 && b == 1.0
        return quadsum(f, points, weights)
    else
        return quadsum(adapt(f, a, b), points, weights)
    end
end