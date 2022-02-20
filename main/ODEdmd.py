from fileinput import filename
import matplotlib.pyplot as plt
import numpy as np

from pydmd import DMD


def check_float(potential_float):
    try:
        float(potential_float)
        return True
    except ValueError:
        return False


def read(file_path):

    file = open(file_path, "r")
    leitura = []
    for line in file:
        lido = line.strip().split(";")
        if check_float(lido[0]) == False:
            leitura.append(lido)
        else:
            leitura.append(lido)

    file.close()

    return leitura


def vectorize(file_path):
    leitura = read(file_path)
    time = []
    pos = []
    for i in range(0, len(leitura)):
        time.append(float(leitura[i][0]))
        pos.append(float(leitura[i][1]))
    return time, pos


t, x = vectorize("data1.txt")


def makeplot(x, t):
    plt.plot(t, x)
    plt.title(
        r"Solution for $\ddot{\theta} + [\delta + \varepsilon\cos{(t)}]\sin{(\theta)} = 0$"
    )
    plt.savefig("images/plot.png")


makeplot(x, t)

xgrid, tgrid = np.meshgrid(x, t)

dmd = DMD(svd_rank=2)
dmd.fit(xgrid.T)

for eig in dmd.eigs:
    print(
        "Eigenvalue {}: distance from unit circle {}".format(
            eig, np.abs(np.sqrt(eig.imag ** 2 + eig.real ** 2) - 1)
        )
    )

dmd.plot_eigs(show_axes=True, show_unit_circle=True, filename="images/eigs.pdf")

print (dmd.eigs)
# for mode in dmd.modes.T:
#    print(mode.real)
#    plt.plot(x, mode.real)
#    plt.title("Modes")
# plt.savefig("images/modes.png")

# for dynamic in dmd.dynamics:
#    plt.plot(t, dynamic.real)
#    plt.title("Dynamics")
# plt.savefig("images/dynamics.png")

# plt.plot(tgrid, xgrid)

# plt.pcolor(xgrid, tgrid)
# plt.colorbar()
# plt.show()

# for mode in dmd.modes.T:
#   plt.plot(x, mode.real)
#  plt.title("Modes")
# plt.savefig("images/modes.png")


# for dynamic in dmd.dynamics:
#    plt.plot(t, dynamic.real)
#    plt.title("Dynamics")
# plt.savefig("images/dynamics.png")
