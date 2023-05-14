import math
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg


def VAR_SECOND(x):
    return math.cos(pow(x, 2))


def get_nodes_values(nodes_count, A, B, function):
    X = [A]
    Y = [function(A)]
    for i in range(nodes_count):
        x = float(input(f"Введите значение x узла {i + 1}: "))
        X.append(x)
        Y.append(function(x))
    X.append(B)
    Y.append(function(B))
    return X, Y


def Lagrange(x, X, Y):
    temp = []
    for i in range(len(X)):
        temp.append(Y[i] * l_Lagrange(i, x, X))
    return sum(temp)


def l_Lagrange(i, x, X):
    numerator = 1
    denominator = 1
    for j in range(len(X)):
        if j != i:
            numerator *= x - X[j]
    for j in range(len(X)):
        if j != i:
            denominator *= X[i] - X[j]
    return float(numerator/denominator)


def spline_solution(x, X, Y):
    h = []
    for i in range(1, len(X)):
        h.append(X[i] - X[i - 1])

    n = len(X) - 1
    koef = []
    for i in range(0, n):
        temp = [0] * 3 * n
        temp[i] = h[i]
        temp[n + i] = pow(h[i], 2)
        temp[2 * n + i] = pow(h[i], 3)
        koef.append(temp)

    for i in range(0, n - 1):
        temp = [0] * 3 * n
        temp[i + 1] = 1
        temp[i] = -1
        temp[n + i] = -2 * h[i]
        temp[2 * n + i] = -3 * pow(h[i], 2)
        koef.append(temp)

    for i in range(0, n - 1):
        temp = [0] * 3 * n
        temp[n + i + 1] = 1
        temp[n + i] = -1
        temp[2 * n + i] = -3 * pow(h[i], 2)
        koef.append(temp)

    temp = [0] * 3 * n
    temp[n] = 1
    koef.append(temp)

    temp = [0] * 3 * n
    temp[2*n - 1] = 1
    temp[3*n - 1] = 3 * h[n - 1]
    koef.append(temp)

    values = [0] * 3 * n
    for i in range(1, n + 1):
        values[i-1] = Y[i] - Y[i - 1]

    n_koef = np.array(koef)
    n_values = np.array(values)
    solve = numpy.linalg.solve(n_koef, n_values)
    a_koefs = [Y[i] for i in range(n)]
    b_koefs = solve[:n]
    c_koefs = solve[n:2*n]
    d_koefs = solve[2*n:3*n]

    temp_x = x
    i = 0
    while temp_x - h[i] > 0:
        temp_x -= h[i]
        i += 1

    result = [a_koefs[i], b_koefs[i] * (x - X[i]), c_koefs[i] * pow(x - X[i], 2), d_koefs[i] * pow(x - X[i], 3)]
    return sum(result)


if __name__ == '__main__':
    COUNT_NODES = 4
    A = 0
    B = 2
    X, Y = get_nodes_values(COUNT_NODES, A, B, VAR_SECOND)

    x = []
    y1 = []
    y2 = []
    y3 = []
    STEPS = 100
    STEP = (B - A) / STEPS
    for i in range(STEPS):
        xi = A + STEP * i
        x.append(xi)
        y1.append(VAR_SECOND(xi))
        y2.append(Lagrange(xi, X, Y))
        y3.append(spline_solution(xi, X, Y))
    plt.plot(x, y1, color="b")
    plt.plot(x, y2, color="g")
    plt.plot(x, y3, color="r")
    plt.show()