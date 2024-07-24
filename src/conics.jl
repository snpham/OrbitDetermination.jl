module Conics

using LinearAlgebra
include("constants.jl")
using .Constants

mutable struct OrbitalElements
    mu::Float64
    rvec::Vector{Float64}
    vvec::Vector{Float64}
    r_mag::Float64
    v_mag::Float64
    r_hat::Vector{Float64}
    v_hat::Vector{Float64}
    h_vec::Vector{Float64}
    h_mag::Float64
    h_hat::Vector{Float64}
    node_vec::Vector{Float64}
    node_mag::Float64
    e_vec::Vector{Float64}
    e_mag::Float64
    e_hat::Vector{Float64}
    oe::Vector{Float64}
end

function OrbitalElements(rvec::Vector{Float64}, vvec::Vector{Float64}, center::String="earth")
    mu = Constants.mu(center)
    
    r_mag = norm(rvec)
    v_mag = norm(vvec)
    r_hat = rvec / r_mag
    v_hat = vvec / v_mag

    h_vec = cross(rvec, vvec)
    h_mag = norm(h_vec)
    h_hat = h_vec / h_mag

    node_vec = cross([0.0, 0.0, 1.0], h_vec)
    node_mag = norm(node_vec)

    e_vec = eccentricity_vector(rvec, vvec, mu, r_mag, v_mag)
    e_mag = norm(e_vec)
    e_hat = e_vec / e_mag

    oe = [
        semimajor_axis(mu, h_mag, e_mag),
        e_mag,
        rad2deg(inclination(h_vec, h_mag)),
        rad2deg(raan(node_vec, node_mag, inclination(h_vec, h_mag))),
        rad2deg(aop(node_vec, node_mag, e_vec, e_mag)),
        rad2deg(true_anomaly(e_vec, e_mag, rvec, r_mag, vvec))
    ]

    return OrbitalElements(mu, rvec, vvec, r_mag, v_mag, r_hat, v_hat, h_vec, h_mag, h_hat, node_vec, node_mag, e_vec, e_mag, e_hat, oe)
end

function eccentricity_vector(rvec, vvec, mu, r_mag, v_mag)
    scalar1 = v_mag^2 / mu - 1.0 / r_mag
    term1 = scalar1 * rvec
    term2 = -dot(rvec, vvec) / mu * vvec
    return term1 + term2
end

function inclination(h_vec, h_mag)
    return acos(h_vec[3] / h_mag)
end

function true_anomaly(e_vec, e_mag, rvec, r_mag, vvec)
    if e_mag == 0
        return NaN
    end
    ta = acos(dot(e_vec, rvec) / (e_mag * r_mag))
    if dot(rvec, vvec) < 0
        ta = 2 * π - ta
    end
    return ta
end

function raan(node_vec, node_mag, inclination)
    if inclination == 0
        return NaN
    end
    omega = acos(node_vec[1] / node_mag)
    if node_vec[2] < 0
        omega = 2 * π - omega
    end
    return omega
end

function aop(node_vec, node_mag, e_vec, e_mag)
    if e_mag == 0 || node_mag == 0
        return NaN
    end
    argp = acos(dot(node_vec, e_vec) / (node_mag * e_mag))
    if e_vec[3] < 0
        argp = 2 * π - argp
    end
    return argp
end

function semiparameter(mu, h_mag)
    return h_mag^2 / mu
end

function semimajor_axis(mu, h_mag, e_mag)
    return semiparameter(mu, h_mag) / (1 - e_mag^2)
end

function energy(v_mag, mu, r_mag)
    return v_mag^2 / 2 - mu / r_mag
end

function period(mu, semimajor_axis)
    return 2 * π * sqrt(semimajor_axis^3 / mu)
end

function rp(semimajor_axis, e_mag)
    return semimajor_axis * (1 - e_mag)
end

function ra(semimajor_axis, e_mag)
    return semimajor_axis * (1 + e_mag)
end

# Example usage (assuming Constants.mu is defined)
# rvec = [7000.0, 0.0, 0.0]
# vvec = [0.0, 7.5, 1.0]
# k = OrbitalElements(rvec, vvec)
# println(k.oe)
end