class Satellite:
    def __init__(self, x, y, z, m, h, v):
        # m
        self.x = x
        # m
        self.y = y
        # m
        self.z = z
        # kg
        self.mass = m
        # m
        self.height = h
        # km/s
        self.velocity = v
        # m2
        self.area = x*z
        # unitless
        self.dragCoeff = 0.8
"""
    def __init__(self, a, e, i, omega, w, theta):
        self.a = a
        self.e = e
        self.i = i
        self.omega = omega
        self.w = w
        self.theta = theta
    
    def __init__(self, x, y, z, dx, dy, dz):
        self.x = x
        self.y = y
        self.z = z
        self.dx = dx
        self.dy = dy
        self.dz = dz
"""
