# test plotting
import matplotlib.pyplot as plt
import numpy as np


def plotFreq():
    lines = open("../out/ctp2.txt").readlines()
    imag = []
    for i in range(len(lines)):
        split = lines[i].split(", ")
        imag.append(split[1])

    npimag = np.array(imag, dtype=float)

    # plt.plot(2 / 8000 * np.abs(npimag[:8000 // 2]))
    plt.plot(np.abs(npimag[:8000 // 2]))

    plt.show()


if __name__ == '__main__':
    plotFreq()
