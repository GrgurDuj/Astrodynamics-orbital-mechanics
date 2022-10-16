import math
from math import sqrt
import numpy as np

import matplotlib.pyplot as plt

import Atmos
from Atmos import airDensity
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
# in kg
earthMass = 5.97219 * 10 ** 24
# in km
earthMRadius = 6378
# in m
h = 500000
sat2Velocity = sqrt(earthMu / (h / 1000 + earthMRadius)) + 0.001
sat2Height = earthMu / (sat2Velocity ** 2)

sat1 = Satellite(0.300, 0.300, 0.300, 30, h, sqrt(earthMu / (h / 1000 + earthMRadius)))
sat2 = Satellite(0.300, 0.300, 0.300, 30, (sat2Height - earthMRadius) * 1000, sat2Velocity)

print(sat2Height - earthMRadius)


def solve_kepler(M, e):
    '''Given some mean anomaly, M (rad), and eccentricity, e, this function
    finds the eccentric anomaly E from the relation M = E - e*sin(E).

    Parameters
    ----------
    M : float
        Mean anomaly (radians) between +/- pi
    e : float
        Osculating eccentricity between (not inclusive) 0 and 1

    Returns
    -------
    E : float
        Eccentric anomaly (radians)

    '''

    E1 = M  # Initialise eccentric anomaly
    residual = 1.0  # Initialise convergence residual
    while residual >= 0.00001:
        fn = E1 - (e * math.sin(E1)) - M
        fd = 1 - (e * math.cos(E1))
        E2 = E1 - (fn / fd)
        residual = abs(E2 - E1)  # Compute residual
        E1 = E2  # Update the eccentric anomaly

    return E2


def orbitPeriod(a):
    return 2 * math.pi * math.sqrt(a ** 3 / earthMu)

def HohmannTransferDV(rd,maintenance_tolerance,maintenance_fro):
    R1 = rd
    if maintenance_fro == True:
        R2 = rd + (2*maintenance_tolerance*1000)
    else:
        R2 = rd + (maintenance_tolerance*1000)
    DV1 = math.sqrt(earthMu/R1) * (math.sqrt((2*R2)/(R1+R2)) - 1)
    DV2 = math.sqrt(earthMu/R2) * (1 - math.sqrt((2*R1)/(R1+R2)))
    totalDV = DV1+DV2
    Time_Elapsed = math.pi * math.sqrt(((R1+R2)**3)/(8*earthMu))
    # Returns (Delta_V (m/s), Time_Elapsed (seconds))
    return totalDV, Time_Elapsed

def orbitalLinearVelocity(h):
    # r is height
    # a is sma
    return sqrt(earthMu / (h / 1000 + earthMRadius))


def dragAcceleration(r, a, sat):
    assert isinstance(sat, Satellite)

    altitude = r - earthMRadius * 1000  # Altitude (m)
    a2m = sat.area / sat.mass  # Area to mass ratio
    velocity = orbitalLinearVelocity(altitude)  # Velocity magnitude (m/s)
    atmDensity = Atmos.airDensity(altitude)  # Density (kg/m^3)
    return 0.5 * atmDensity * sat.dragCoeff * a2m * (velocity ** 2)


def altitudeDecayRate(r, a, sat):
    dragAccel = dragAcceleration(r, a, sat)
    return -1 * (dragAccel * orbitPeriod(a) / math.pi)


def calculateDistance(anomaly1, anomaly2, altitude1, altitude2):
    anomalyDifference = anomaly2 - anomaly1
    tendonLength = 2 * (altitude2 + earthMRadius) * np.sin(anomalyDifference / 2)
    altitudeDifference = altitude1 - altitude2
    distance = sqrt(altitudeDifference ** 2 + tendonLength ** 2)  # √[(x₂ - x₁)² + (y₂ - y₁)²]

    return distance


def keplerToCartesian(a, e, i, RAAN, omega, nu, mu):
    p = a * (1 - e ** 2)  # semi-latus rectum, [km]
    r = p / (1 + e * np.cos(nu))  # orbit radius, [km]
    # h = np.sqrt(mu*a*(1-e^2)) # angular momentum
    x = r * (np.cos(RAAN) * np.cos(omega + nu) - np.sin(RAAN) * np.sin(omega + nu) * np.cos(i))  # x-position, [km]
    y = r * (np.sin(RAAN) * np.cos(omega + nu) + np.cos(RAAN) * np.sin(omega + nu) * np.cos(i))  # y-position, [km]
    z = r * (np.sin(i) * np.sin(omega + nu))  # z-position, [km]
    cart = [x, y, z]  # cartesian state vector
    return cart


def calculate(sat):
    orb_altitude = sat.height
    orb_eccentricity = 0
    orb_mean_anomaly = 0
    orb_semimajor = sat.height + earthMRadius * 1000
    orb_inclination = 90
    orb_true_anomaly = 0
    orb_w = 0
    orb_omega = 0
    altitudes = []
    mean_anomalies = []

    """
    have to keep track of
    shape size 
    eccentricity e          ✓ doesn't change with decay
    semi-major axis a       ✓

    orbital plane
    inclination i           ✓ doesn't change with decay
    longitude of ascending node omega   ✓ doesn't change with decay

    argument of periapsis w - ascending node to periapsis   ✓ doesn't change with decay
    mean anomaly v - position at time ✓ 
    """
    decay_alt = 0.0  # The decay amount in meters
    # total_DV = 0.0 # Total Delta-V used (m/s)
    # thrustcount = 0 # Counting the number of thrusts needed
    meanAnomaly = np.deg2rad(orb_mean_anomaly)  # Init the mean anomaly in radians
    # cart_state = keplerToCartesian(orb_a, orb_e, orb_i, orb_rightA, orb_w, orb_m, earthMu)
    lifetime_flag = False
    total_duration = 0  # seconds
    tstep = 1
    current_altitude = orb_altitude
    while lifetime_flag == False and total_duration < 60 * 60 * 24 * 7:
        # update time step
        total_duration = total_duration + tstep

        # mean motion
        keplerPeriod = orbitPeriod(current_altitude + earthMRadius)
        meanMotion = (2 * np.pi) / keplerPeriod

        # mean anomaly at time step
        meanAnomaly = (meanAnomaly + np.pi + (meanMotion * tstep))
        meanAnomaly = meanAnomaly % (2 * math.pi) - math.pi
        mean_anomalies.append(meanAnomaly)

        # eccentric anomaly
        eccnAnomaly = solve_kepler(meanAnomaly, orb_eccentricity)

        # Now, we can substitute orb_E to find the satellite radial distance.
        rd = orb_semimajor * (1 - orb_eccentricity * np.cos(eccnAnomaly))

        # compute decay rate and apply to semimajor axis
        decay_rate = altitudeDecayRate(rd, orb_semimajor, sat)
        decay_alt = decay_rate * tstep
        orb_semimajor += decay_alt

        current_altitude = rd / 1000 - earthMRadius
        altitudes.append(current_altitude)

        if current_altitude < 100.0:
            lifetime_flag = True
            lifetime = total_duration

    if lifetime_flag:
        lifetime_days = total_duration / 86400
        lifetime_orbits = total_duration / orbitPeriod(sat.height)
        lifetime_str = "Lifetime decay is"
        lifetime_str += 'after ' + str(lifetime_orbits) + ' orbits.'
        lifetime_str += 'The lifetime is ' + str(lifetime_days) + ' days. \n'
        print(lifetime_str)

    if not lifetime_flag:
        print("Satellite did not decay in 180 days")

    lifetime_orbits = total_duration / orbitPeriod(sat.height)
    # print("Lifetime orbits: " + str(lifetime_orbits))
    print(total_duration / 60 / 60 / 24)
    return altitudes, mean_anomalies


major_ticks = np.arange(0, 60 * 24 * 366, 60 * 24 * 30)
major_labels = np.arange(0, 13, 1)

print(str(sat1.velocity) + "km/s")
print(str(sat2.velocity) + "km/s")
altitudesSat2, meanAnomaliesSat2 = calculate(sat2)
"""
plt.plot(altitudesSat2)
print(altitudesSat2[-1])
print(sat2.height / 1000 - altitudesSat2[-1])
plt.ylim(ymin=498.18, ymax=498.195)
plt.title("Satellite 2 decay")
plt.xticks(major_ticks, labels=major_labels)
plt.xlabel("Days since orbit start")
plt.savefig("Sat2decay.png")
plt.show()
"""
altitudesSat1, meanAnomaliesSat1 = calculate(sat1)
"""
print(altitudesSat1[-1])
print(sat1.height / 1000 - altitudesSat1[-1])
plt.plot(altitudesSat1)
plt.ylim(ymin=499.99, ymax=500.001)
plt.title("Satellite 1 decay")
plt.xticks(major_ticks, labels=major_labels)
plt.xlabel("Days since orbit start")
plt.savefig("Sat1decay.png")
# plt.xticks(minor_ticks, labels=None, minor=True)
plt.show()

plt.plot(meanAnomaliesSat2)
plt.title("Mean anomalies 2")
plt.xticks(major_ticks, labels=major_labels)
plt.xlabel("Days since orbit start")
plt.savefig("meanAnom2total.png")
plt.show()

plt.plot(meanAnomaliesSat1)
plt.title("Mean anomalies 1")
plt.xticks(major_ticks, labels=major_labels)
plt.xlabel("Days since orbit start")
plt.savefig("meanAnom1total.png")
plt.show()

plt.plot(meanAnomaliesSat2)
plt.xlim(xmin=0, xmax=2 * 60 * 60)
plt.xticks(np.arange(0, 3 * 60 * 60, 60 * 60), np.arange(0, 3, 1))
plt.xlabel("Hours since orbit start")
plt.title("mean anomaly sat 2")
plt.savefig("meanAnom22hours.png")
plt.show()

plt.plot(meanAnomaliesSat1)
plt.xlim(xmin=0, xmax=2 * 60 * 60)
plt.title("mean anomaly sat 1")
plt.xticks(np.arange(0, 3 * 60 * 60, 60 * 60), np.arange(0, 3, 1))
plt.xlabel("Hours since orbit start")
plt.savefig("meanAnom12hours.png")
plt.show()
"""
def calcVelocity(height):
    return sqrt(earthMu / (height + earthMRadius))


orbitPeriods1 = []
orbitPeriods2 = []
velocities1 = []
velocities2 = []

for alt in altitudesSat1:
    velocities1.append(calcVelocity(alt))
    orbitPeriods1.append(orbitPeriod(alt))

for alt in altitudesSat2:
    velocities2.append(calcVelocity(alt))
    orbitPeriods2.append(orbitPeriod(alt))


"""
print("Sat 1: min" + str(sat1.velocity) + "sat 1: max " + str(velocities1[-1]) + "diff: " + str(
    sat1.velocity - velocities1[-1]))
plt.plot(velocities1)
plt.title("velocities 1")
plt.xticks(major_ticks, labels=major_labels)
plt.xlabel("Days since orbit start")
plt.ylim(ymin=7.61268, ymax=7.61269)
plt.savefig("velocities1.png")
plt.show()

print("Sat 2: min" + str(sat2.velocity) + "sat 2: max " + str(velocities2[-1]) + "diff: " + str(
    sat2.velocity - velocities2[-1]))
plt.plot(velocities2)
plt.title("velocities 2")
plt.xticks(major_ticks, labels=major_labels)
plt.xlabel("Days since orbit start")
plt.ylim(ymin=7.61368, ymax=7.61369)
plt.savefig("velocities2.png")
plt.show()
"""

differences = []
for (v1,v2) in zip(velocities1, velocities2):
    differences.append(v2-v1)
plt.plot(differences)
plt.xticks(major_ticks, labels=major_labels)
plt.xlabel("Days since orbit start")
plt.title("Velocity differences")
plt.show()

plt.plot(orbitPeriods1)
plt.xticks(major_ticks, labels=major_labels)
plt.xlabel("Days since orbit start")
plt.title("Satellite 1 orbit periods")
plt.show()

plt.plot(orbitPeriods2)
plt.xticks(major_ticks, labels=major_labels)
plt.xlabel("Days since orbit start")
plt.title("Satellite 2 orbit periods")
plt.show()

odifferences = []
for (o1, o2) in zip(orbitPeriods1, orbitPeriods2):
    odifferences.append(o2 - o1)
plt.plot(odifferences)
plt.title("Orbit period difference")
plt.xlabel("Days since orbit start")
plt.title("Satellite orbit period differences")
plt.xticks(major_ticks, labels=major_labels)
plt.show()

"""
distances = []
count = 0
for (an1, an2, alti1, alti2) in zip(meanAnomaliesSat1, meanAnomaliesSat2, altitudesSat1, altitudesSat2):
    count = + 1
    a1, a2, alt1, alt2 = (an1, an2, alti1, alti2)
    distances.append(calculateDistance(a1, a2, alt1, alt2))
np.savetxt("distances.csv", distances)
plt.plot(distances)
plt.title("Distance")
plt.xticks(major_ticks, labels=major_labels)
plt.xlabel("Months   since orbit start")
plt.savefig("distances.png")
plt.show()
"""

print("Sat 1 init: " + str(sat1.velocity) + " Sat 1 last: " + str(velocities1[-1]) + " Sat 1 change: " + str(sat1.velocity - velocities1[-1]))
print("Sat 2 init: " + str(sat2.velocity) + " Sat 2 last: " + str(velocities2[-1]) + " Sat 2 change: " + str(sat2.velocity - velocities2[-1]))