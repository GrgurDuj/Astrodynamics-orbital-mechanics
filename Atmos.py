from scipy.interpolate import interp1d

def airDensity(h):
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
    heights = [0, 100, 150, 200, 250, 300, 350, 400, 450, 500]

    interpolated_densities = interp1d(heights, densities)
    h = h / 1000
    return interpolated_densities(h)
