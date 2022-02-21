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

def makeplot(x, t,name):
    plt.plot(t, x)
    plt.title("Solution for $\ddot{\theta} + [\delta + \varepsilon\cos{(t)}]\sin{(\theta)} = 0$")
    plt.savefig(name)
    print("\n" + name + " has been generated\n")

def main():
    t, x = vectorize("data1.txt")

    makeplot(x,t,"plot.png")

    xgrid, tgrid = np.meshgrid(x, t)

    dmd = DMD(svd_rank=2)
    dmd.fit(xgrid.T)

    print ("\nEigenValues: " + str(dmd.eigs))

    for eig in dmd.eigs:
        print("\nEigenvalue " + str(eig) + ": distance from unit circle " +  str(np.abs(np.sqrt(eig.imag ** 2 + eig.real ** 2) - 1)))
        print("Parte Real: " + str(eig.real) + " Parte Imaginaria: " + str(eig.imag))

    dmd.plot_eigs(show_axes=True, show_unit_circle=True, filename="images/eigs.png")

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

main()