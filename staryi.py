import matplotlib.pyplot as plt
import numpy as np
from scipy.special import comb
from functools import reduce
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
# import antonutils.geometry as geom
import time
import warnings

warnings.simplefilter('ignore', np.RankWarning)

dotnum = 100
t_otn = 0.6
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111)  # , projection='3d', adjustable='box')
ax.set_autoscaley_on(False)


def tochka_peresecheniya(verhnaya_liniya, nizhnaya_liniya):
    verhnaya_liniyax, verhnaya_liniyay = verhnaya_liniya
    nizhnaya_liniyax2, nizhnaya_liniyay2 = nizhnaya_liniya
    A_1 = verhnaya_liniyay[0] - verhnaya_liniyay[1]
    B_1 = verhnaya_liniyax[1] - verhnaya_liniyax[0]
    C_1 = verhnaya_liniyax[0] * verhnaya_liniyay[1] - \
          verhnaya_liniyax[1] * verhnaya_liniyay[0]
    A_2 = nizhnaya_liniyay2[0] - nizhnaya_liniyay2[1]
    B_2 = nizhnaya_liniyax2[1] - nizhnaya_liniyax2[0]
    C_2 = nizhnaya_liniyax2[0] * nizhnaya_liniyay2[1] - \
          nizhnaya_liniyax2[1] * nizhnaya_liniyay2[0]
    if A_1 * B_2 - A_2 * B_1 == 0:
        print('Ошибка. Прямые не пересекаются')
        tochka_peresecheniyax = (verhnaya_liniyax[1] + nizhnaya_liniyax2[0]) / 2
        tochka_peresecheniyay = (verhnaya_liniyay[1] + nizhnaya_liniyay2[0]) / 2
        return [tochka_peresecheniyax, tochka_peresecheniyay]
    else:
        tochka_peresecheniyax = -(C_1 * B_2 - C_2 * B_1) / (A_1 * B_2 - A_2 * B_1)
        tochka_peresecheniyay = -(A_1 * C_2 - A_2 * C_1) / (A_1 * B_2 - A_2 * B_1)
        return [tochka_peresecheniyax, tochka_peresecheniyay]


def bernstein_poly(i, n, t):
    s = comb(n, i) * (t ** (n - i)) * (1 - t) ** i
    return s


def bezier_curve(_points_, ntimes=dotnum):
    n_points = len(_points_)
    x_points = np.array([p[0] for p in _points_])
    z_points = np.array([p[1] for p in _points_])

    t = np.linspace(0.0, 1.0, ntimes)

    polynomial_array = np.array([bernstein_poly(i, n_points - 1, t) for i in range(0, n_points)])

    _xvals = np.dot(x_points, polynomial_array)
    _yvals = np.dot(z_points, polynomial_array)

    return _xvals, _yvals


do_plot = False
plt.style.use('mpl15')
t = np.linspace(0, 2 * np.pi, 200)

beta_10 = 57.418364162592106
beta_11 = 29.19834792437215
beta_1sr = 39.47345303684428
beta_20 = 90.18011833453723
beta_21 = 90.18011833453723
beta_2sr = 90.18011833453723

# в миллиметрах
d_10 = 0.04764912383531226 * 1000
d_11 = 0.13341754673887432 * 1000
d_1sr = 0.09053333528709329 * 1000
# d_20 =  0.19059649534124903*1000
# d_21 =  0.19059649534124903*1000
# d_2sr =  0.19059649534124903*1000
d_20 = d_10 * 1.05
d_21 = d_11 * 0.95
d_2sr = d_1sr

# b_2 =  0.011874967899489663*1000
b_2 = 0

w_10 = 178.0151232783465
w_1sr = 235.95270004400393
w_11 = 307.4813351537137
w_20 = 159.4450173766263
w_2sr = 159.4450173766263
w_21 = 159.4450173766263

# Угол радиального наклона потока на входе в колесо
degree_10 = 0.
degree_11 = 0.
degree_1sr = 0.

# Угол радиального наклона потока на выходе из колесо
# degree_20 = 90.
# degree_21 = 90.
# degree_2sr = 90.
degree_20 = 2.
degree_21 = -2.
degree_2sr = 0.
# Угол раскрытия лопатки
gamma_1 = 2.  # град
# Радиус входной окружности
radius_1 = 2.
# Коээфициент удлинннеия кромки
coeff_1 = 1.
# Угол закрытия лопатки
gamma_2 = 2.  # град
# Радиус входной окружности
radius_2 = 2.
# Коээфициент удлинннеия кромки
coeff_2 = 1.

# Координаты точки начала координат

x_nach_coord, y_nach_coord, z_nach_coord = 0., 0., 0.
smeshnie_x, smeshnie_y, smeshnie_z = 0., 0., 0.
# Длина рабочего колеса по shroud
dlina = d_11

X0 = x_nach_coord + smeshnie_x
Y0 = y_nach_coord + smeshnie_y
Z0 = z_nach_coord + smeshnie_z


def create_interp_lopatka_sechenie(x0, y0, w_1, beta_1, w_2, beta_2, dlina, b_2, gamma_1, radius_1, coeff_1, gamma_2,
                                   radius_2, coeff_2):
    # Координаты точки рабочего колеса на входе
    x_1 = x0
    y_1 = y0
    # Координаты нижней точки рабочего колеса на входе с учетом скорости и угла наклона потока
    x_1w = x0 - w_1 * np.cos(np.radians(beta_1))
    y_1w = y0 - w_1 * np.sin(np.radians(beta_1))
    # Для окружности
    x_1w_radius = x0 - radius_1 * np.cos(np.radians(beta_1)) * coeff_1 * 1.75
    y_1w_radius = y0 - radius_1 * np.sin(np.radians(beta_1)) * coeff_1 * 1.75
    # Координаты нижней точки рабочего колеса на входе
    x_1_spinka = x0 + radius_1 * np.sin(np.radians(beta_1))
    y_1_spinka = y0 - radius_1 * np.cos(np.radians(beta_1))
    # Координаты нижней точки рабочего колеса на входе с учетом скорости и угла наклона потока
    x_1w_spinka = x_1_spinka - w_1 * np.cos(np.radians(beta_1 - gamma_1 / 2))
    y_1w_spinka = y_1_spinka - w_1 * np.sin(np.radians(beta_1 - gamma_1 / 2))
    # Для окружности
    x_1w_spinka_radius = x_1_spinka - radius_1 * np.cos(np.radians(beta_1 - gamma_1 / 2)) * coeff_1
    y_1w_spinka_radius = y_1_spinka - radius_1 * np.sin(np.radians(beta_1 - gamma_1 / 2)) * coeff_1
    # Координаты нижней точки рабочего колеса на входе
    x_1_korytse = x0 - radius_1 * np.sin(np.radians(beta_1))
    y_1_korytse = y0 + radius_1 * np.cos(np.radians(beta_1))
    # Координаты нижней точки рабочего колеса на входе с учетом скорости и угла наклона потока
    x_1w_korytse = x_1_korytse - w_1 * np.cos(np.radians(beta_1 + gamma_1 / 2))
    y_1w_korytse = y_1_korytse - w_1 * np.sin(np.radians(beta_1 + gamma_1 / 2))
    # Для окружности
    x_1w_korytse_radius = x_1_korytse - radius_1 * np.cos(np.radians(beta_1 + gamma_1 / 2)) * coeff_1
    y_1w_korytse_radius = y_1_korytse - radius_1 * np.sin(np.radians(beta_1 + gamma_1 / 2)) * coeff_1

    # Установочные углы
    beta_y = (beta_1 + beta_2) / 2
    # Координаты нижней точки рабочего колеса на входе
    x_2 = x0 + dlina
    y_2 = y0 + dlina * np.tan(np.radians(beta_y))
    # Координаты нижней точки рабочего колеса на входе с учетом скорости и угла наклона потока
    x_2w = x0 + dlina + w_2 * np.cos(np.radians(beta_2))
    y_2w = y0 + dlina * np.tan(np.radians(beta_y)) + w_2 * np.sin(np.radians(beta_2))
    # Для окружности
    x_2w_radius = x0 + dlina + radius_2 * np.cos(np.radians(beta_2)) * coeff_2 * 1.75
    y_2w_radius = y0 + dlina * np.tan(np.radians(beta_y)) + radius_2 * np.sin(np.radians(beta_2)) * coeff_2 * 1.75
    # Координаты нижней точки рабочего колеса на входе
    x_2_spinka = x0 + dlina + radius_2 * np.sin(np.radians(beta_2 + gamma_2 / 2))
    y_2_spinka = y0 + dlina * np.tan(np.radians(beta_y)) - radius_2 * np.cos(np.radians(beta_2 + gamma_2 / 2))
    # Координаты нижней точки рабочего колеса на входе с учетом скорости и угла наклона потока
    x_2w_spinka = x_2_spinka + w_2 * np.cos(np.radians(beta_2 + gamma_2 / 2))
    y_2w_spinka = y_2_spinka + w_2 * np.sin(np.radians(beta_2 + gamma_2 / 2))
    # Для окружности
    x_2w_spinka_radius = x_2_spinka + radius_2 * np.cos(np.radians(beta_2 + gamma_2 / 2)) * coeff_2
    y_2w_spinka_radius = y_2_spinka + radius_2 * np.sin(np.radians(beta_2 + gamma_2 / 2)) * coeff_2

    # Координаты нижней точки рабочего колеса на входе
    x_2_korytse = x0 + dlina - radius_2 * np.sin(np.radians(beta_2 + gamma_2 / 2))
    y_2_korytse = y0 + dlina * np.tan(np.radians(beta_y)) + radius_2 * np.cos(np.radians(beta_2 + gamma_2 / 2))
    # Координаты нижней точки рабочего колеса на входе с учетом скорости и угла наклона потока
    x_2w_korytse = x_2_korytse + w_2 * np.cos(np.radians(beta_2 - gamma_2 / 2))
    y_2w_korytse = y_2_korytse + w_2 * np.sin(np.radians(beta_2 - gamma_2 / 2))
    # Для окружности
    x_2w_korytse_radius = x_2_korytse + radius_2 * np.cos(np.radians(beta_2 - gamma_2 / 2)) * coeff_2
    y_2w_korytse_radius = y_2_korytse + radius_2 * np.sin(np.radians(beta_2 - gamma_2 / 2)) * coeff_2

    tchk_1 = tochka_peresecheniya([[x_1w, x_1], [y_1w, y_1]], [[x_2w, x_2], [y_2w, y_2]])
    tchk_1_spinka = tochka_peresecheniya(
        [[x_1w_spinka, x_1_spinka],
         [y_1w_spinka, y_1_spinka]],
        [[x_2w_spinka, x_2_spinka],
         [y_2w_spinka, y_2_spinka]])
    tchk_1_korytse = tochka_peresecheniya(
        [[x_1w_korytse, x_1_korytse],
         [y_1w_korytse, y_1_korytse]],
        [[x_2w_korytse, x_2_korytse],
         [y_2w_korytse, y_2_korytse]])

    points = [
        [x_1, y_1],
        [tchk_1[0], tchk_1[1]],
        [x_2, y_2]
    ]
    points_spinka = [
        [x_1_spinka, y_1_spinka],
        [tchk_1_spinka[0], tchk_1_spinka[1]],
        [x_2_spinka, y_2_spinka]
    ]
    points_korytse = [
        [x_1_korytse, y_1_korytse],
        [tchk_1_korytse[0], tchk_1_korytse[1]],
        [x_2_korytse, y_2_korytse]
    ]
    points_1radius = [
        [x_1_korytse, y_1_korytse],
        [x_1w_korytse_radius, y_1w_korytse_radius],
        [x_1w_radius, y_1w_radius],
        [x_1w_spinka_radius, y_1w_spinka_radius],
        [x_1_spinka, y_1_spinka]
    ]
    points_2radius = [
        [x_2_korytse, y_2_korytse],
        [x_2w_korytse_radius, y_2w_korytse_radius],
        [x_2w_radius, y_2w_radius],
        [x_2w_spinka_radius, y_2w_spinka_radius],
        [x_2_spinka, y_2_spinka]
    ]
    xvals, yvals = bezier_curve(points)
    xvals_spinka, yvals_spinka = bezier_curve(points_spinka)
    xvals_korytse, yvals_korytse = bezier_curve(points_korytse)
    xvals_1radius, yvals_1radius = bezier_curve(points_1radius, ntimes=20)
    xvals_2radius, yvals_2radius = bezier_curve(points_2radius, ntimes=20)

    if do_plot:
        fig = plt.figure(figsize=(16, 9))
        plt.plot(x_1, y_1, 'ro')
        plt.plot(x_1w, y_1w, 'ro')
        plt.plot([x_1w, x_1], [y_1w, y_1], '--')
        plt.plot(x_1w_radius, y_1w_radius, 'ro')
        plt.plot([x_1w_radius, x_1], [y_1w_radius, y_1], '--')

        plt.plot(x_2, y_2, 'o', color='c')
        plt.plot(x_2w, y_2w, 'o', color='c')
        plt.plot([x_2w, x_2], [y_2w, y_2], '--')
        plt.plot(x_2w_radius, y_2w_radius, 'o', color='c')
        plt.plot([x_2w_radius, x_2], [y_2w_radius, y_2], '--')

        plt.plot(x_1_spinka, y_1_spinka, 'ro')
        plt.plot(x_1w_spinka, y_1w_spinka, 'ro')
        plt.plot([x_1w_spinka, x_1_spinka], [y_1w_spinka, y_1_spinka], '--', color='grey')
        plt.plot(x_1w_spinka_radius, y_1w_spinka_radius, 'ro')
        plt.plot([x_1w_spinka_radius, x_1_spinka], [y_1w_spinka_radius, y_1_spinka], '--', color='grey')

        plt.plot(x_2_spinka, y_2_spinka, 'o', color='c')
        plt.plot(x_2w_spinka, y_2w_spinka, 'o', color='c')
        plt.plot([x_2w_spinka, x_2_spinka], [y_2w_spinka, y_2_spinka], '--', color='grey')
        plt.plot(x_2w_spinka_radius, y_2w_spinka_radius, 'o', color='c')
        plt.plot([x_2w_spinka_radius, x_2_spinka], [y_2w_spinka_radius, y_2_spinka], '--', color='grey')

        plt.plot(x_1_korytse, y_1_korytse, 'ro')
        plt.plot(x_1w_korytse, y_1w_korytse, 'ro')
        plt.plot([x_1w_korytse, x_1_korytse], [y_1w_korytse, y_1_korytse], '--', color='grey')
        plt.plot(x_1w_korytse_radius, y_1w_korytse_radius, 'ro')
        plt.plot([x_1w_korytse_radius, x_1_korytse], [y_1w_korytse_radius, y_1_korytse], '--', color='grey')

        plt.plot(x_2_korytse, y_2_korytse, 'o', color='c')
        plt.plot(x_2w_korytse, y_2w_korytse, 'o', color='c')
        plt.plot([x_2w_korytse, x_2_korytse], [y_2w_korytse, y_2_korytse], '--', color='grey')
        plt.plot(x_2w_korytse_radius, y_2w_korytse_radius, 'o', color='c')
        plt.plot([x_2w_korytse_radius, x_2_korytse], [y_2w_korytse_radius, y_2_korytse], '--', color='grey')

        plt.plot(tchk_1[0], tchk_1[1])
        plt.plot([tchk_1[0]], [tchk_1[1]], 'o')

        plt.plot(tchk_1_spinka[0], tchk_1_spinka[1])
        plt.plot([tchk_1_spinka[0]], [tchk_1_spinka[1]], 'o')

        plt.plot(tchk_1_korytse[0], tchk_1_korytse[1])
        plt.plot([tchk_1_korytse[0]], [tchk_1_korytse[1]], 'o')

        plt.plot(xvals, yvals)
        plt.plot(xvals_spinka, yvals_spinka)
        plt.plot(xvals_korytse, yvals_korytse)
        plt.plot(xvals_1radius, yvals_1radius)
        plt.plot(xvals_2radius, yvals_2radius)

        plt.axis('equal')
        plt.grid()
        plt.show()

    coords = [
        np.concatenate((np.flip(xvals_2radius), xvals_spinka, xvals_1radius, np.flip(xvals_korytse))),
        # xvals_2radius+xvals_spinka+xvals_1radius+xvals_korytse,
        np.concatenate((np.flip(yvals_2radius), yvals_spinka, yvals_1radius, np.flip(yvals_korytse))),
        # yvals_2radius+yvals_spinka+yvals_1radius+yvals_korytse
    ]
    # f1 = np.ones(len(coords[0]))
    # f2 = np.ones(len(coords[1]))
    center_of_mass = [np.sum(coords[0] * np.ones(len(coords[0]))) / np.sum(np.ones(len(coords[0]))),
                      np.sum(coords[1] * np.ones(len(coords[1]))) / np.sum(np.ones(len(coords[1])))]
    # plt.plot(coords[0], coords[1])
    # plt.plot(center_of_mass[0], center_of_mass[1], 'ro')
    profil_x = coords[0] - center_of_mass[0]
    profil_y = coords[1] - center_of_mass[1]

    return [profil_x, profil_y]


profil_x1, profil_y1 = create_interp_lopatka_sechenie(
    X0, X0, w_10, beta_10, w_20, beta_20, dlina, b_2,
    gamma_1, radius_1, coeff_1, gamma_2, radius_2, coeff_2
)

# print(profil_x1)
# print(profil_y1)

plt.plot(profil_x1, profil_y1, 'k')


# plt.plot(profil_x1, profil_y1)
# plt.show()

index_of_max_y = np.argmax(profil_y1)
index_of_min_y = np.argmin(profil_y1)

max_x = profil_x1[index_of_max_y]
min_x = profil_x1[index_of_min_y]
max_y = profil_y1[index_of_max_y]
min_y = profil_y1[index_of_min_y]

delta_y = max_y - min_y

profil_x2, profil_y2 = profil_x1 + (t_otn * (max_y - min_y)), profil_y1

plt.plot(profil_x2, profil_y2, 'k')
# plt.plot([max_x], [max_y], 'ro')
# plt.plot([min_x], [min_y], 'ro')
plt.plot([-delta_y / 2, delta_y], [max_y, max_y], color='gray')
plt.plot([-delta_y / 2, delta_y], [min_y, min_y], color='gray')

ggggg = (len(profil_x1) == len(profil_y1))

# Спинка первого профиля
splitx_sp = profil_x1[np.argmax(profil_y1):np.argmin(profil_y1)]
splity_sp = profil_y1[np.argmax(profil_y1):np.argmin(profil_y1)]

# Корытце второй профиля
splitx_kor = np.delete(profil_x2, range(np.argmax(profil_y1) - 0, np.argmin(profil_y1) - 0))
splity_kor = np.delete(profil_y2, range(np.argmax(profil_y1) - 0, np.argmin(profil_y1) - 0))
# print(np.argmin(profil_y), np.argmax(profil_y))
# print(len(splitx_sp), len(splity_sp), len(splitx_kor), len(splity_kor), len(profil_x), len(profil_y))

splitx_kor_fixed = [i for i, y in sorted(zip(splitx_kor, splity_kor), key=lambda pair: pair[1])]
splity_kor_fixed = [y for i, y in sorted(zip(splitx_kor, splity_kor), key=lambda pair: pair[1])]


def diff(m, n):
    return np.sqrt((n[0] - m[0]) ** 2 + (n[1] - m[1]) ** 2)


# Средняя линия тока
plt.plot([max_x], [max_y], 'ro')
plt.plot([min_x], [min_y], 'ro')
y_range = np.linspace(max(min(splity_sp), min(splity_kor_fixed)), min(max(splity_sp), max(splity_kor_fixed)), 1000)

interpol1 = interpolate.interp1d(splity_sp, splitx_sp)  # np.interp(y_range, splity_sp, splitx_sp)
interpol2 = interpolate.interp1d(splity_kor_fixed, splitx_kor_fixed)
# np.interp(y_range, splity_kor_fixed, splitx_kor_fixed)

# plt.plot(interpol1(y_range), y_range, 'c')
# plt.plot(interpol2(y_range), y_range, 'c')

XXXXXX = (interpol2(y_range) + interpol1(y_range)) / 2
YYYYYY = y_range
XXXXXX2 = (interpol2(y_range) + interpol1(y_range)) / 2 - np.mean(interpol2(y_range) + interpol1(y_range))
YYYYYY2 = y_range

plt.plot((interpol2(y_range) + interpol1(y_range)) / 2, y_range, 'c')
# plt.axis('equal')
# plt.show()

# y_circl = np.linspace(min(max(splity_sp), max(splity_kor_fixed)), max(min(splity_sp), min(splity_kor_fixed)), 10)
# x_circl = (interpol2(y_circl) + interpol1(y_circl))/2

# for xx, yy in zip(x_circl, y_circl):
#     plt.plot([xx], [yy], 'ko')

b = 1
dlina_spinka_list = []
dlina_kor_list = []


def length(x: list, y: list):
    ret = []
    for i in range(len(x)):
        if i < len(x) - 1:
            ret.append(
                np.sqrt((x[i + 1] - x[i]) ** 2 + (y[i + 1] - y[i]) ** 2)
            )
    return sum(ret)


dlina_spinka = length(splitx_sp, splity_sp)
dlina_kor = length(splitx_kor_fixed, splity_kor_fixed)

percentage = [i * 10 for i in np.linspace(0.0, 9, 30)]
dlina_spinka_per = [u / 100 * dlina_spinka for u in percentage]
dlina_kor_per = [u / 100 * dlina_kor for u in percentage]

import time

start = time.time()
i_ret = []

for ppp in dlina_spinka_per:
    dlina_spinka_list = []
    for i in range(len(splitx_sp)):
        # bernstein_i += int(i_ret/2.2)  # +i_ret*ppp*dlina_spinka*0.98
        # print(bernstein_i)
        if i % 100 == 0:
            # print(bernstein_i)
            pass
        if i < len(splitx_sp) - 1:
            dlina_spinka_list.append(
                np.sqrt((splitx_sp[i + 1] - splitx_sp[i]) ** 2 + (splity_sp[i + 1] - splity_sp[i]) ** 2)
            )
        else:
            pass
        # print(sum(dlina_spinka_list))
        if sum(dlina_spinka_list) >= ppp:
            # print('OK')
            # plt.plot([splitx_sp[bernstein_i]], [splity_sp[bernstein_i]], 'ko')
            i_ret.append(i)
            break
        else:
            pass
# print('Время расчета', time.time() - start)
# normals = []
# print('i_ret', i_ret)

list_of_cirl_linesx1 = []
list_of_cirl_linesy1 = []
list_of_cirl_linesx2 = []
list_of_cirl_linesy2 = []

for num in i_ret:
    if num < len(splitx_sp) - 1:
        list_x = [splitx_sp[num - 1], splitx_sp[num], splitx_sp[num + 1]]
        list_y = [splity_sp[num - 1], splity_sp[num], splity_sp[num + 1]]
    elif num == 0:
        list_x = [splitx_sp[num], splitx_sp[num + 1], splitx_sp[num + 2]]
        list_y = [splity_sp[num], splity_sp[num + 1], splity_sp[num + 2]]
    else:
        list_x = [splitx_sp[num - 1], splitx_sp[num], splitx_sp[0]]
        list_y = [splity_sp[num - 1], splity_sp[num], splity_sp[0]]

    polinom = np.poly1d(np.polyfit(list_x, list_y, 2))
    # print('polinom', polinom)
    deriv_polinom = polinom.deriv()
    # print('deriv_polinom', deriv_polinom)
    k_x = deriv_polinom(splitx_sp[num])
    # print('k_x', k_x)
    alpha_k = np.arctan(k_x)
    y = lambda x: -1 / k_x * (x - splitx_sp[num]) + polinom(splitx_sp[num])
    xxx = lambda y: k_x * (polinom(splitx_sp[num]) - y) + splitx_sp[num]
    '''
    # Построение окружностей
    prev = -1.
    for xx, yy in zip(x_circl, y_circl):
        if yy == y_circl[2]:
            # print('y(xx)', y(xx), 'yy', yy,
            #       'np.around(y(xx), 2) == np.around(yy, 2)', (yy-y(xx))/yy, int(y(xx)) == int(yy))
            pass
        current = (yy-y(xx))/yy
        # print('np.sign(current) == np.sign(prev)', np.sign(current) == np.sign(prev))
        # print('np.sign(current)', np.sign(current))
        # print('np.sign(prev)', np.sign(prev))
        if np.sign(current) == prev:
            prev = np.sign(current)
            pass
        else:
            print('np.sign(current) == np.sign(prev)', np.sign(current) == np.sign(prev))
            print('np.sign(current)', np.sign(current))
            print('np.sign(prev)', np.sign(prev))
            print('num', num)
            rangex = np.linspace(splitx_sp[num], splitx_sp[num]+70, 100)
            # plt.plot(rangex, y(rangex), ':')
            plt.plot([xx], [y(xx)], 'ro')
            prev = -1.
    '''
    # print(num, i_ret, num in i_ret)
    rangex = np.linspace(splitx_sp[num], splitx_sp[num] + 70, 100)
    # plt.plot(rangex, y(rangex), ':')
    # normals.append([rangex, y(rangex)])
    # print('diff', range_x - profil_x)

    y_range = np.linspace(
        max(min(splity_sp),
            min(splity_kor_fixed)),
        min(max(splity_sp),
            max(splity_kor_fixed)),
        20000)

    first = (interpol2(y_range) + interpol1(y_range)) / 2
    # plt.plot(first, y_range, 'k')
    second = k_x * (polinom(splitx_sp[num]) - y_range) + splitx_sp[num]
    # plt.plot(second, y_range, 'k--')
    # plt.plot(rangex, y(rangex), 'k--')
    dot_num_sp = np.argwhere(np.diff(np.sign(second - first))).flatten()
    # dot_num_kor = np.argwhere(np.diff(np.sign(kor_y_interp(linspace_kory) - x_kasat_kor))).flatten()
    # print('количество точек пересечения', len(dot_num_sp))
    for i in dot_num_sp:
        # plt.plot(second[bernstein_i], y_range[bernstein_i], 'bo')
        # for bernstein_i in dot_num_kor:
        #     plt.plot(kor_y_interp(linspace_kory)[bernstein_i], linspace_kory[bernstein_i], 'bo')
        if min(splity_sp) <= y_range[i] <= max(splity_sp):
            # Построение окружности:
            # plt.plot(splitx_sp[num], splity_sp[num], 'ko')
            rad = length([splitx_sp[num], second[i]], [splity_sp[num], y_range[i]])
            # rad2 = length([splitx_sp[num], second[bernstein_i]], [splity_sp[num], y_range[bernstein_i]])
            if abs((rad-(max_y - min_y))/rad) < 1:
                pass
            else:
                # print('rad', rad)
                # print((rad-(max_y - min_y))/rad)
                plt.plot(rad * np.sin(t) + second[dot_num_sp[0]],
                         rad * np.cos(t) + y_range[dot_num_sp[0]], ':',
                         color='gray')
                plt.plot(rad * np.sin(t) + (t_otn * (max_y - min_y))*2.5,
                         rad * np.cos(t) + y_range[dot_num_sp[0]],
                         'k:')
                list_of_cirl_linesx1.append(rad * np.sin(np.pi/2) + (t_otn * (max_y - min_y))*2.5)
                list_of_cirl_linesy1.append(rad * np.cos(np.pi/2) + y_range[dot_num_sp[0]])
                list_of_cirl_linesx2.append(rad * np.sin(-np.pi/2) + (t_otn * (max_y - min_y))*2.5)
                list_of_cirl_linesy2.append(rad * np.cos(-np.pi/2) + y_range[dot_num_sp[0]])

plt.plot(list_of_cirl_linesx1, list_of_cirl_linesy1, 'c')
plt.plot(list_of_cirl_linesx2, list_of_cirl_linesy2, 'c')

i_ret2 = []
for ppp in dlina_kor_per:
    dlina_spinka_list = []
    for i in range(len(splitx_kor_fixed)):
        # bernstein_i += int(i_ret/2.2)  # +i_ret*ppp*dlina_spinka*0.98
        # print(bernstein_i)
        if i % 100 == 0:
            # print(bernstein_i)
            pass
        if i < len(splitx_kor_fixed) - 1:
            dlina_spinka_list.append(
                np.sqrt((splitx_kor_fixed[i + 1] - splitx_kor_fixed[i]) ** 2 +
                        (splity_kor_fixed[i + 1] - splity_kor_fixed[i]) ** 2)
            )
        else:
            pass

        if sum(dlina_spinka_list) >= ppp:
            # plt.plot([splitx_kor_fixed[bernstein_i]], [splity_kor_fixed[bernstein_i]], 'ko')
            # plt.text(splitx_kor_fixed[bernstein_i], splity_kor_fixed[bernstein_i], str(bernstein_i))
            i_ret2.append(i)
            break
        else:
            pass
# print('Время расчета', time.time() - start)
# normals = []
for num in i_ret2:
    if num < len(splitx_kor_fixed) - 1:
        list_x2 = [splitx_kor_fixed[num - 1], splitx_kor_fixed[num], splitx_kor_fixed[num + 1]]
        list_y2 = [splity_kor_fixed[num - 1], splity_kor_fixed[num], splity_kor_fixed[num + 1]]
    elif num == 0:
        list_x2 = [splitx_kor_fixed[num], splitx_kor_fixed[num + 1], splitx_kor_fixed[num + 2]]
        list_y2 = [splity_kor_fixed[num], splity_kor_fixed[num + 1], splity_kor_fixed[num + 2]]
    else:
        list_x2 = [splitx_kor_fixed[num - 1], splitx_kor_fixed[num], splitx_kor_fixed[0]]
        list_y2 = [splity_kor_fixed[num - 1], splity_kor_fixed[num], splity_kor_fixed[0]]

    polinom2 = np.poly1d(np.polyfit(list_x2, list_y2, 2))
    # print('polinom', polinom)
    deriv_polinom2 = polinom2.deriv()
    # print('deriv_polinom', deriv_polinom)
    k_x2 = deriv_polinom2(splitx_kor_fixed[num])
    # print('k_x2', k_x2)
    alpha_k2 = np.arctan(k_x2)
    y2 = lambda x: -1 / k_x2 * (x - splitx_kor_fixed[num]) + polinom2(splitx_kor_fixed[num])
    # Построение окружностей
    # for xx, yy in zip(x_circl, y_circl):
    #     if np.around(y2(xx), 4) == np.around(yy, 4):
    #         rangex2 = np.linspace(splitx_kor_fixed[num] - 70, splitx_kor_fixed[num], 100)
    #         plt.plot(rangex2, y2(rangex2), ':')

    if num in i_ret2:
        # print(num, i_ret2, num in i_ret2)
        # normals.append([rangex, y(rangex)])
        # plt.plot(rangex2, y2(rangex2), ':')
        pass

# print(normals)
plt.axis('equal')
plt.xlim(-150, (t_otn * (max_y - min_y))*4)
plt.ylim(-150, 150)
plt.grid()
plt.show()
del(fig)
