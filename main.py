from math import sqrt

from Satellite import Satellite

earthmu = 398600.440
earthmradius = 6371
sat1 = Satellite(0.300, 0.300, 0.300, 30, 500000, sqrt(earthmu / 500000))
print(sat1.v)