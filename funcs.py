
from scipy.special import comb
import numpy as np


def tochka_peresecheniya(verhnaya_liniyax, verhnaya_liniyay, nizhnaya_liniyax2, nizhnaya_liniyay2):
    a_1 = verhnaya_liniyay[0] - verhnaya_liniyay[1]
    b_1 = verhnaya_liniyax[1] - verhnaya_liniyax[0]
    c_1 = verhnaya_liniyax[0] * verhnaya_liniyay[1] - verhnaya_liniyax[1] * verhnaya_liniyay[0]
    a_2 = nizhnaya_liniyay2[0] - nizhnaya_liniyay2[1]
    b_2 = nizhnaya_liniyax2[1] - nizhnaya_liniyax2[0]
    c_2 = nizhnaya_liniyax2[0] * nizhnaya_liniyay2[1] - nizhnaya_liniyax2[1] * nizhnaya_liniyay2[0]
    if a_1 * b_2 - a_2 * b_1 == 0:
        # print('Ошибка. Прямые не пересекаются')
        tochka_peresecheniyax = (verhnaya_liniyax[1] + nizhnaya_liniyax2[0]) / 2
        tochka_peresecheniyay = (verhnaya_liniyay[1] + nizhnaya_liniyay2[0]) / 2
        return [tochka_peresecheniyax, tochka_peresecheniyay]
    else:
        tochka_peresecheniyax = -(c_1 * b_2 - c_2 * b_1) / (a_1 * b_2 - a_2 * b_1)
        tochka_peresecheniyay = -(a_1 * c_2 - a_2 * c_1) / (a_1 * b_2 - a_2 * b_1)
        return [tochka_peresecheniyax, tochka_peresecheniyay]


def bernstein_poly(i, n, t):
    s = comb(n, i) * (t ** (n - i)) * (1 - t) ** i
    return s


def bezier_curve(_points_, ntimes=200):
    n_points = len(_points_)
    x_points = np.array([p[0] for p in _points_])
    z_points = np.array([p[1] for p in _points_])

    t = np.linspace(0.0, 1.0, ntimes)

    polynomial_array = np.array([bernstein_poly(i, n_points - 1, t) for i in range(0, n_points)])

    _xvals = np.dot(x_points, polynomial_array)
    _yvals = np.dot(z_points, polynomial_array)

    return _xvals, _yvals
