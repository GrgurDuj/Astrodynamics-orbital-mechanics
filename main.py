from math import sqrt
from math import exp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from Satellite import Satellite

"""
• Spacecraft body size is X = 300 mm, Y= 300 mm & Z=300 mm.
• Spacecraft total mass is 30 kg.
• The spacecraft is pointing nadir with its Z-axis, and in the flight direction with its X-axis.
• The orbit is sun-synchronous, circular at 500 km altitude.
• The solar array is body mounted.
• The satellites are ejected in the same orbit with one of the two satellites having +1 m/s
velocity in the flight direction w.r.t the other satellite which is exactly in the reference orbit.
• Each satellite is oriented with the X-Z plane perpendicular to the flight direction.
• The minimum thrust level of the propulsion unit is 1 N, the maximum thrust level is 1 kN.
"""
# in km s
earthMu = 398600.440
# in km
earthMRadius = 6371
# in m
h = 500000
heights = [0, 100, 150, 200, 250, 300, 350, 400, 450, 500]
densities = [1.225,
             5.25 * 10 ** -7,
             1.73 * 10 ** -9,
             2.41 * 10 ** -10,
             5.97 * 10 ** -11,
             1.87 * 10 ** -11,
             6.66 * 10 ** -12,
             2.62 * 10 ** -12,
             1.09 * 10 ** -12,
             4.76 * 10 ** -13]
interpolated_densities = interp1d(heights, densities)

sat1 = Satellite(0.300, 0.300, 0.300, 30, h, sqrt(earthMu / (h / 1000 + earthMRadius)))
sat2 = Satellite(0.300, 0.300, 0.300, 30, h, sqrt(earthMu / (h / 1000 + earthMRadius)) + 0.001)
print(str(sat1.v) + "km/s")
print(str(sat2.v) + "km/s")

"""
f10 min is 70
f10 max is 300
T = 900 + 2.5 (F10 – 70) + 1.5 Ap	(Kelvin)
μ = 27 – 0.012 (h – 200)	180 < h (km) < 500
H = T / μ	(km)
ρ = 6x10-10 exp( - (h –175) / H )	( kg m-3)
"""
def airDensity(interpolated_densities, h):
    h = h / 1000
    return interpolated_densities(h)
    """h = h / 1000
    f10 = 300
    temperature = 900 + 2.5 * (f10 - 70)
    mu = 27 - 0.012 * (h - 200)
    H = temperature / mu
    density = (6 * 10 ** -10) ** (-(h - 175) / H)
    # temperature = -131.21 + 0.00299 * h
    # density = 2.488 * ((temperature + 273.1)/216.6) ** -11.388
    return density
    """


print(airDensity(interpolated_densities, h))
"""
Calculates air density at specified eight in meters
"""
