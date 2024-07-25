module Conics

using LinearAlgebra
include("constants.jl")
include("rotations.jl")
include("matrices.jl")
using .Constants
using .Rotations
using .Matrix
using Base.Math: sin, cos, sqrt

mutable struct OrbitalElements
    sma::Float64
    e_mag::Float64
    incl::Float64
    raan::Float64
    argp::Float64
    ta::Float64
    r_hat::Vector{Float64}
    v_hat::Vector{Float64}
    h_hat::Vector{Float64}
    e_hat::Vector{Float64}
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

    sma = semimajor_axis(mu, h_mag, e_mag)
    incl = rad2deg(inclination(h_vec, h_mag))
    rasc = rad2deg(raan(node_vec, node_mag, inclination(h_vec, h_mag)))
    argp = rad2deg(aop(node_vec, node_mag, e_vec, e_mag))
    ta = rad2deg(true_anomaly(e_vec, e_mag, rvec, r_mag, vvec))

    return OrbitalElements(sma, e_mag, incl, rasc, argp, ta, r_hat, v_hat, h_hat, e_hat)
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

function get_rv_frm_elements(elements, center::String="earth", method::String="sma")
    """
    Computes position/velocity vectors from orbital elements.
    We first compute pos/vel in the PQW system, then rotate to the
    geocentric equatorial system.

    Parameters:
    elements - for method = 'p': (p, e, i, raan, aop, ta)
               for method = 'sma': (a, e, i, raan, aop, ta)
    a: semi-major axis (km); e: eccentricity
    i: inclination (rad); raan: right ascending node (rad)
    aop: argument of periapsis (rad); ta: true anomaly (rad)
    center: center object of orbit; default=earth

    Returns:
    rvec: positional vectors of spacecraft [IJK] (km)
    vvec: velocity vectors of spacecraft [IJK] (km/s)
    """

    if method == "p"
        p, e, i, raan, aop, ta = elements
    elseif method == "sma"
        a, e, i, raan, aop, ta = elements
        p = a * (1 - e^2)
    else
        throw(ArgumentError("Invalid method. Use 'p' or 'sma'."))
    end

    # Determine which planet center to compute from
    mu = Constants.mu(center)

    if method == "sma"
        r = p / (1 + e * cos(ta))
        h = sqrt(mu * a * (1 - e^2))
        rvec = [
            r * (cos(raan) * cos(aop + ta) - sin(raan) * sin(aop + ta) * cos(i)),
            r * (sin(raan) * cos(aop + ta) + cos(raan) * sin(aop + ta) * cos(i)),
            r * (sin(i) * sin(aop + ta))
        ]
        vvec = [
            rvec[1] * h * e * sin(ta) / (r * p) - h / r * (cos(raan) * sin(aop + ta) + sin(raan) * cos(aop + ta) * cos(i)),
            rvec[2] * h * e * sin(ta) / (r * p) - h / r * (sin(raan) * sin(aop + ta) - cos(raan) * cos(aop + ta) * cos(i)),
            rvec[3] * h * e * sin(ta) / (r * p) + h / r * (sin(i) * cos(aop + ta))
        ]
        return vcat(rvec, vvec)

    elseif method == "p"
        # Assigning temporary variables
        aop_t = aop
        raan_t = raan
        ta_t = ta 

        # Checking for undefined states 
        if e == 0 && i == 0
            aop_t = 0.0
            raan_t = 0.0
            ta_t = aop_t + raan_t + ta
        elseif e == 0
            aop_t = 0.0
            ta_t = aop_t + ta
        elseif i == 0
            raan_t = 0.0
            aop_t = raan_t + aop
            ta_t = ta
        end

        # Converting elements into state vectors in PQW frame
        r_pqw = [p * cos(ta_t) / (1 + e * cos(ta_t)), p * sin(ta_t) / (1 + e * cos(ta_t)), 0.0]
        v_pqw = [-sqrt(mu / p) * sin(ta_t), sqrt(mu / p) * (e + cos(ta_t)), 0.0]
        
        # Get 313 transformation matrix to geocentric-equatorial frame
        m1 = Rotations.rotate(-aop, 'z')
        m2 = Rotations.rotate(-i, 'x')
        m3 = Rotations.rotate(-raan, 'z')
        T_ijk_pqw = Matrix.mxm(m3, Matrix.mxm(m2, m1))

        # state vector from PQW to ECI
        r_ijk = Matrix.mxv(T_ijk_pqw, r_pqw)
        v_ijk = Matrix.mxv(T_ijk_pqw, v_pqw)
        return vcat(r_ijk, v_ijk)
    else
        throw(ArgumentError("Invalid method. Use 'p' or 'sma'."))
    end
end

# Example usage (assuming Constants.mu is defined)
# rvec = [7000.0, 0.0, 0.0]
# vvec = [0.0, 7.5, 1.0]
# k = OrbitalElements(rvec, vvec)
# println(k.oe)
end