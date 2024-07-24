module Constants

# gravitational coefficients, km3/s2
const planetary_mu = Dict(
    "sun" => 1.32712440018e11,
    "mercury" => 2.2032e13,
    "venus" => 3.24858599e5,
    "earth" => 3.986004415e5,
    "moon" => 4.9028000661e3,
    "mars" => 4.28283100e4,
    "jupiter" => 1.266865361e8,
    "saturn" => 3.7931208e7,
    "uranus" => 5.7939513e6,
    "neptune" => 6.835100e6,
    "pluto" => 8.71e2,
)

const mu_sun = planetary_mu["sun"]
const mu_venus = planetary_mu["venus"]
const mu_earth = planetary_mu["earth"]
const mu_moon = planetary_mu["moon"]
const mu_mars = planetary_mu["mars"]
const mu_jupiter = planetary_mu["jupiter"]

# universal gravitational constant
const G_univ = 6.67408e-20  # km3/(kg)

# periods, days
const planetary_period = Dict(
    "mercury" => missing,
    "venus" => 224.701,
    "earth" => 365.242189,
    "moon" => missing,
    "mars" => missing,
    "jupiter" => missing,
    "saturn" => missing,
    "uranus" => missing,
    "neptune" => missing,
    "pluto" => missing,
)
const P_earth = 365.242189  # days
const P_venus = 224.701  # days
const AU = 1.49597870691e8  # km

# planetary radii, km
const planetary_radii = Dict(
    "mercury" => missing,
    "venus" => 6051.8,
    "earth" => 6378.14,
    "moon" => 1737.4,
    "mars" => 3396.19,
    "jupiter" => 71492,
    "saturn" => 60268,
    "uranus" => 25559,
    "neptune" => 24764,
    "pluto" => 1188.3,
)
const r_venus = 6051.8
const r_earth = 6378.14
const r_moon = 1737.14
const r_mars = 3396.19
const r_jupiter = 71492

# planetary heliocentric semi-major axes (km)
const planetary_sma = Dict(
    "mercury" => 57.909e6,
    "venus" => 0.723 * AU,  # from wikipedia
    "earth" => 149598023,
    "moon" => 384400,  # from earth
    "mars" => 227939186,  # 1.523679342 AU
    # from nssdc.gsfc.nasa.gov/planetary/factsheet/
    "jupiter" => 778.570e6,
    "saturn" => 1433.529e6,
    "uranus" => 2872.463e6,
    "neptune" => 4495.060e6,
    "pluto" => 5869.656e6,
)

const sma_earth = 149598023
const sma_mars = 227939186  # 1.523679342 AU
const sma_jupiter = 778.570e6  # from nssdc.gsfc.nasa.gov/planetary/factsheet/
const sma_moon = 384399

# earth constants
const E_earth = 0.081819221456
const omega_earth = 7.292115e-5  # +/- 1.5e-12 (rad/s)
const mass_earth = 5.9722e24  # km

function mu(center::String="earth")
    planet = lowercase(center)
    if haskey(planetary_mu, planet)
        return planetary_mu[planet]
    else
        throw(ArgumentError("Undefined planet, ending.."))
    end
end

function sma(center::String="earth")
    planet = lowercase(center)
    if haskey(planetary_sma, planet)
        return planetary_sma[planet]
    else
        throw(ArgumentError("Undefined planet, ending.."))
    end
end

function period(center::String="earth")
    planet = lowercase(center)
    if haskey(planetary_period, planet)
        return planetary_period[planet]
    else
        throw(ArgumentError("Undefined planet, ending.."))
    end
end

function radius(center::String="earth")
    planet = lowercase(center)
    if haskey(planetary_radii, planet)
        return planetary_radii[planet]
    else
        throw(ArgumentError("Undefined planet, ending.."))
    end
end

end # module
