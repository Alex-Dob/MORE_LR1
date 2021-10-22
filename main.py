import numpy as np
import scipy.stats as sps


def FindNormRaspr(Z):
    print('\nРазбиение массива на классы')
    k = 1 + 3.32 * np.log10(n)
    delta = Z.max() - Z.min()
    step = delta / np.round(k)
    CLASS = np.linspace(Z.min(), Z.max(), num=int(np.round(k)))
    print('k = ', np.round(k))
    print('Размах варьирования = ', delta)

    X = np.zeros(int(np.round(k)) - 1)
    B = np.zeros(int(np.round(k)) - 1)

    for i in range(len(CLASS) - 1):
        for j in range(n):
            if Z[j] >= CLASS[i] and Z[j] <= CLASS[i + 1]:
                B[i] = B[i] + 1
        avg = (CLASS[i] + CLASS[i + 1]) / 2
        X[i] = avg

    print('\n\t\tИнтервалы', '\t\t Середины интервалов', '\t Частота')
    for i in range(len(CLASS) - 1):
        print('[', np.round(CLASS[i], 5), ':', np.round(CLASS[i + 1], 5), ']\t',
              np.round(X[i], 5), '\t\t\t\t\t', int(B[i]))

    average = sum(B * X) / n
    s_ = np.sqrt((sum(B * X ** 2) - sum(B * X) ** 2 / n) / (n - 1))
    v = s_ / average
    k_ = n * step / s_
    z = (X - average) / s_

    print('\naverage = ', average)
    print('s_ = ', s_)
    print('v = ', v)
    print('k_ = ', k_)
    print('z = ', z)

    E = k_ * (1 / (2 * np.pi) * np.exp(-z ** 2 / 2))

    XI_2 = sum((B - E) ** 2 / E)

    print('\nE = ', E)

    xi2 = sps.chi2.ppf(1 - 0.1, int(k) - 1 - 2)

    if XI_2 < xi2:
        print('\nГипотеза нормальности распределения может быть принята на 10 %-ном уровне, т. к. XI_2 < xi2')
    else:
        print('\nГипотеза нормальности распределения отвергается, т. к. XI_2 > xi2')

    F_E = np.zeros(int(np.round(k)) - 1)
    F_B = np.zeros(int(np.round(k)) - 1)

    for i in range(len(F_E)):
        if i == 0:
            F_E[i] = E[i]
            F_B[i] = B[i]
        else:
            F_E[i] = E[i] + F_E[i - 1]
            F_B[i] = B[i] + F_B[i - 1]

    D = np.max(F_B - F_E) / n
    D_ = 0.218
    if D < D_:
        print('Гипотеза нормальности распределения может быть принята на 10 %-ном уровне, т. к. D < D_')
    else:
        print('Гипотеза нормальности распределения отвергается, т. к. D > D_')


Z = np.genfromtxt('Variant_04.txt', delimiter="\n")
DEL = []
print(Z)
n = 0
s2 = _s = s = average = 0
m1 = m3 = m4 = 0
v = 0
t1 = t2 = 0
tau = 1
while t1 < tau:
    n = len(Z)
    average = Z.mean()
    s2 = Z.var()
    s = np.sqrt(s2)
    _s = np.sqrt(sum((Z - average) ** 2) / (n - 1))
    m1 = sum(Z - average) / n
    m3 = sum((Z - average) ** 3) / n
    m4 = sum((Z - average) ** 4) / n
    v = _s / average
    d_max = abs(Z - average).max()
    tau = d_max / _s

    print("average = ", average)
    print("s2 =", s2)
    print("s =", s)
    print("_s =", _s)
    print("m1 =", m1)
    print("m3 =", m3)
    print("m4 =", m4)
    print("v =", v)
    print("d_max =", d_max)
    print("tau =", tau)

    alpha_1 = 0.05
    alpha_2 = 0.001

    t_1 = sps.t.ppf(1 - alpha_1, n - 2)
    t_2 = sps.t.ppf(1 - alpha_2, n - 2)
    print('t_1 = ', t_1)
    print('t_2 = ', t_2)

    t1 = t_1 * np.sqrt(n - 1) / np.sqrt(n - 2 + t_1 ** 2)
    t2 = t_2 * np.sqrt(n - 1) / np.sqrt(n - 2 + t_1 ** 2)
    print('t1 = ', t1)
    print('t2 = ', t2)

    if t1 < tau or t2 < tau:
        x = abs(Z[0] - average)
        k = 0
        for i in range(n):
            if x < abs(Z[i] - average):
                x = abs(Z[i] - average)
                k = i
        print('Удаляем элемент: ', Z[k])
        DEL = np.append(DEL, Z[k])
        Z = np.delete(Z, k)
    print('\n')


print("Отсеянные элементы:")
print(DEL)
print("\nПосле отсева грубой погршености получим:")
print(np.sort(Z))

g1 = g2 = 0
G1 = G2 = 0
S_G1 = S_G2 = 0
print("\nПроверка нормальности распределения")
if v * 100 < 33:
    g1 = m3 / s2 ** (3 / 2)
    g2 = m4 / s2 ** 2 - 3
    G1 = np.sqrt(n - 1) / (n - 2) * g1
    G2 = (n - 1) * ((n + 1) * g2 + 6) / ((n - 2) * (n - 3))
    S_G1 = np.sqrt(6 * n * (n - 1) / ((n - 2) * (n + 1) * (n + 3)))
    S_G2 = np.sqrt(24 * n * (n - 1) ** 2 / ((n - 3) * (n - 2) * (n + 3) * (n + 5)))

    print('g1 = ', g1)
    print('g2 = ', g2)
    print('G1 = ', G1)
    print('G2 = ', G2)
    print('S_G1 = ', S_G1)
    print('S_G2 = ', S_G2)

if abs(G1) <= 3 * S_G1 and abs(G2) <= 5 * S_G2:
    print('\n', np.round(abs(G1), 2), ' <= ', np.round(3 * S_G1, 2),
          '\tи\t', np.round(abs(G2), 2), ' <= ', np.round(5 * S_G2, 2))
    print('- соблюдение этих условий говорит о возможности ',
          'принятия гипотезы нормального распрераспределения')

#Z = 1 / np.sqrt(Z)
#Z = 1 / Z
#Z = Z ** 1.5
#Z = Z ** 2
Z = np.log10(Z + 0.1) * 10 ** 0.1
FindNormRaspr(Z)
