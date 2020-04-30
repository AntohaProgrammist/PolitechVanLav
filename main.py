import numpy as np

# Хорда
B = 30  # мм

# Угол в относительном движении в сечении перед рабочей лопаткой
beta_1 = 38  # градусов
# Угол в относительном движении в сечении за рабочей лопаткой
beta_2 = 34  # градусов
# Установочный угол
beta_y = 180 - (beta_1+beta_2)

# Угол раскрытия лопатки
gamma_1 = 8  # градусов
# Угол закртия лопатки
gamma_2 = 4  # градусов

# Скорости в относительном движении перед и за рабочим колесом турбины
w_1 = 314/100  # м/с
w_2 = 371/100  # м/с

# Радиусы скргуления входной и выходной кромок
radius_1 = 0.04 * B
radius_2 = 0.01 * B

# Координаты начала построений
X, Y = 0, 0

# Перевод в радианы
beta_1 = np.radians(beta_1)
beta_2 = np.radians(beta_2)
beta_y = np.radians(beta_y)
gamma_1 = np.radians(gamma_1)
gamma_2 = np.radians(gamma_2)

# Координаты точки рабочего колеса на входе
x_1 = X
y_1 = Y
# Координаты н и угла наклона потижней точки рабочего колеса на входе с учетом скоростиока
x_1w = x_1 + w_1 * np.cos(beta_1)
y_1w = y_1 + w_1 * np.sin(beta_1)

# Это точка l
x_1_spinka = x_1 - radius_1 * np.sin(beta_1-gamma_1/2)
y_1_spinka = y_1 + radius_1 * np.cos(beta_1-gamma_1/2)
x_1w_spinka = x_1_spinka + w_1 * np.cos(beta_1-gamma_1/2)
y_1w_spinka = y_1_spinka + w_1 * np.sin(beta_1-gamma_1/2)

# Это точка g
x_1_korytse = x_1 + radius_1 * np.sin(beta_1+gamma_1/2)
y_1_korytse = y_1 - radius_1 * np.cos(beta_1+gamma_1/2)
x_1w_korytse = x_1_korytse + w_1 * np.cos(beta_1+gamma_1/2)
y_1w_korytse = y_1_korytse + w_1 * np.sin(beta_1+gamma_1/2)

x_2 = x_1 - B*np.cos(beta_y)
y_2 = y_1 - B*np.sin(beta_y)

print('Расстоение между входом и выходом. Хорда = ', np.sqrt((x_2-x_1)**2+(y_2-y_1)**2))

x_2w = x_2 + w_2*np.cos(beta_2)
y_2w = y_2 - w_2*np.sin(beta_2)

x_2_spinka = x_2 - radius_2*np.sin(beta_2-gamma_2/2)
y_2_spinka = y_2 - radius_2*np.cos(beta_2-gamma_2/2)

x_2w_spinka = x_2_spinka + w_2*np.cos(beta_2-gamma_2/2)
y_2w_spinka = y_2_spinka - w_2*np.sin(beta_2-gamma_2/2)

x_2_korytse = x_2 + radius_2*np.sin(beta_2+gamma_2/2)
y_2_korytse = y_2 + radius_2*np.cos(beta_2+gamma_2/2)

x_2w_korytse = x_2_korytse + w_2*np.cos(beta_2+gamma_2/2)
y_2w_korytse = y_2_korytse - w_2*np.sin(beta_2+gamma_2/2)

import matplotlib.pyplot as plt

plt.plot([x_1], [y_1], 'o')
plt.plot([x_1, x_1w], [y_1, y_1w], 'k')
plt.plot([x_1_spinka, x_1w_spinka], [y_1_spinka, y_1w_spinka], 'b--')
plt.plot([x_1_korytse, x_1w_korytse], [y_1_korytse, y_1w_korytse], 'r--')

plt.plot([x_2], [y_2], 'o')
plt.plot([x_2, x_2w], [y_2, y_2w], 'k')
plt.plot([x_2_spinka, x_2w_spinka], [y_2_spinka, y_2w_spinka], 'b--')
plt.plot([x_2_korytse, x_2w_korytse], [y_2_korytse, y_2w_korytse], 'r--')

# plt.plot([x_1w_korytse_radius, x_1w_korytse], [y_1w_korytse_radius, y_1w_korytse], label='1 в корытце')
# plt.plot([x_1w_spinka_radius, x_1w_spinka], [y_1w_spinka_radius, y_1w_spinka], label='1 на спинке')

# plt.plot([x_2], [y_2], 'o', label='2 в корытце')
# # plt.plot([x_2w_spinka_radius, x_2w_spinka], [y_2w_spinka_radius, y_2w_spinka], label='2 на спинке')
# plt.plot([x_2w_korytse_radius, x_2w_korytse], [y_2w_korytse_radius, y_2w_korytse], label='2 в корытце')
# plt.plot([x_2w_spinka_radius, x_2w_spinka], [y_2w_spinka_radius, y_2w_spinka], label='2 на спинке')


plt.title('Название')
plt.xlabel('Ось ИКС')
plt.ylabel('Ось Y')
plt.grid()
plt.axis('equal')
# plt.legend()
plt.show()

b = 1
