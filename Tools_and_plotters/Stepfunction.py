import numpy as np
import matplotlib.pyplot as plt


def step(x, x1, x2):
    y = np.zeros(len(x))

    for i in range(0, len(x)):
        if min(x1, x2) < x[i] < max(x1, x2):
            y[i] = 1

    return y


x = np.arange(0, 100, 0.1)
y = step(x, 25, 50)

plt.plot(x, y)
plt.show()