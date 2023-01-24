import numpy as np
from colorama import Fore
import matplotlib.pyplot as plt
from numpy import linalg


def main():
    def start_coordinates(case):
        x = []
        y = []
        dx = p / (number - 1)
        for i in range(number):
            x.append(i * dx)
            y.append(0.0)

        if case == 0:
            return x
        else:
            return y

    def force(x, y, case):
        force_x = np.zeros(len(x))
        force_y = np.zeros(len(y))

        for i in range(len(x) - 1):
            distance = np.sqrt((x[i + 1] - x[i]) ** 2 + (y[i + 1] - y[i]) ** 2)

            if (distance - length_spring) > 0:
                springs_f_x = c * (distance - length_spring) * (x[i + 1] - x[i]) / (length_spring * distance)
                springs_f_y = c * (distance - length_spring) * (y[i + 1] - y[i]) / (length_spring * distance)

                force_x[i] += springs_f_x
                force_y[i] += (springs_f_y - m * g)
                force_x[i + 1] -= springs_f_x
                force_y[i + 1] -= springs_f_y
            else:
                force_x[i] += 0
                force_y[i] += 0
                force_x[i + 1] -= 0
                force_y[i + 1] -= 0

        force_y[len(x) - 1] -= m * g

        if case == 0:
            return force_x
        else:
            return force_y

    def velocity(x, y, vx, vy, case, fstng, b):
        velocity_x = np.zeros(len(x))
        velocity_y = np.zeros(len(y))

        force_x = force(x, y, 0)
        force_y = force(x, y, 1)

        i = 1
        while i < len(x):
            velocity_x[i] = dt * (force_x[i] - b * vx[i]) / m + vx[i]
            velocity_y[i] = dt * (force_y[i] - b * vy[i]) / m + vy[i]
            i = i + 1

        if fstng <= 3:
            velocity_x[i - 1] = 0
            velocity_y[i - 1] = 0

        if case == 0:
            return velocity_x
        else:
            return velocity_y

    def coordinate(x, y, vx, vy, case):
        for i in range(len(x)):
            x[i] += vx[i] * dt
            y[i] += vy[i] * dt

        x[0] = 0
        y[0] = 0

        if case == 0:
            return x
        else:
            return y

    c = 200  # stiffness
    m = 1  # particle mass
    g = 9.81
    dt = 0.01
    length_spring = 0.05  # segment length
    # dt = 0.05 / (c * length_spring / (m * g))  # time step
    t_1 = dt  # start time
    t_2 = 0
    number = 20  # number of splits
    p = 1  # total length
    b = 10  # viscosity
    fastening = 0  # pinning counter
    eps = 10**(-3)

    print(Fore.BLUE + f'time step: {dt}')

    x_0 = start_coordinates(0)
    y_0 = start_coordinates(1)

    x_1 = x_0.copy()
    y_1 = y_0.copy()

    vx_1 = np.zeros(len(x_1))
    vy_1 = np.zeros(len(y_1))

    while x_1[number-1] > 0:

        yb = 0.0
        xb = - 0.5

        vx_2 = velocity(x_1, y_1, vx_1, vy_1, 0, fastening, b)
        vy_2 = velocity(x_1, y_1, vx_1, vy_1, 1, fastening, b)

        x_2 = coordinate(x_1, y_1, vx_2, vy_2, 0)
        y_2 = coordinate(x_1, y_1, vx_2, vy_2, 1)

        vx_m = linalg.norm(vx_2, np.inf)
        vy_m = linalg.norm(vy_2, np.inf)

        if vx_m < eps and vy_m < eps:
            fastening += 1
            if fastening >= 3:
                b = 0

        if vx_2[number-1] != 0:
            t_2 += dt
            yb = - (1/2) * g * (t_2)**2

        x_1 = x_2.copy()
        y_1 = y_2.copy()

        vx_1 = vx_2.copy()
        vy_1 = vy_2.copy()

        # plotting graph

        plt.axes(xlim=(-p, 1.5 * p), ylim=(-length_spring * number * number / 8, 0.5))
        plt.grid(True)

        plt.plot(x_0, y_0, 'g*-')
        plt.plot(x_1, y_1, 'r*-')
        plt.plot(xb, yb, 'b*')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title(f"Start position and t = {round(t_1, 3)}")

        plt.draw()
        plt.gcf().canvas.flush_events()

        t_1 += dt

        if x_1[number-1] > 0:
            plt.clf()

    print(Fore.CYAN + f'ball_y = {-(1 / 2) * g * (t_2) ** 2}')
    print(Fore.CYAN + f'point_y = {y_2[number - 1]}')


plt.ion()
main()
plt.ioff()
plt.show()
