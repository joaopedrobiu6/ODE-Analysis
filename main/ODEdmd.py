from fileinput import filename
import matplotlib.pyplot as plt
import numpy as np

from pydmd import DMD
from pydmd import HODMD

import pandas as pd  

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
    plt.title(r"Solution for $\ddot{\theta} + [\delta + \varepsilon\cos{(t)}]\sin{(\theta)} = 0$")
    plt.savefig(name)
    print("\n" + name + " has been generated\n")

def dmd_reconstruction(dataframe):
    print(dataframe)
    # create and fit a HODMD model
    hodmd = HODMD(svd_rank=0,
            opt=True,
            exact=True,
            rescale_mode=None,
            forward_backward=False,
            sorted_eigs='abs',
            tlsq_rank=0,d=700)

    hodmd.fit(dataframe.values)

    print("\nEIGENVALUES: (Ordenados por Magnitude)\n")
    print(hodmd.eigs)

    # plot predictions vs actuals on training dataset
    plt.plot(dataframe.values,'.', label = 'Data' , color = 'black', markersize = 3)
    plt.plot(hodmd.reconstructed_data[0].real,label = 'Reconstruction', color = "red")
    plt.legend()
    plt.title("DMD Reconstruction", fontsize=12)
    plt.xlabel(r"${t}$", fontsize=12)
    plt.ylabel(r"${\theta}$", fontsize=12)
    plt.savefig("images/reconstruction.png")

    print("\nFoi criado um plot (images/reconstruction.png) com a reconstrução do sinal pelo DMD\n")

def main():
    #t, x = vectorize("data1.txt")

    #makeplot(x,t,"images/plot.png")

    #Baseado Nisto https://github.com/IvanPMorenoMarcos/DMD-HODM-for-UTSF/blob/main/hodmd.ipynb
    dataframe = pd.read_csv('../ODE-Analysis/data1.txt', index_col=0, header=None, squeeze=True)

    dataframe.plot()
    plt.savefig("images/dataplot.png")
    plt.clf()

    #Quando a solução é estavel a reconstrução funciona bem, quando não é estável mais ou menos
    dmd_reconstruction(dataframe)

    '''
    xgrid, tgrid = np.meshgrid(x, t)
    dmd = DMD(svd_rank=2)
    dmd.fit(xgrid.T)
    print ("\nEigenValues: " + str(dmd.eigs))
    for eig in dmd.eigs:
        print("\nEigenvalue " + str(eig) + ": distance from unit circle " +  str(np.abs(np.sqrt(eig.imag ** 2 + eig.real ** 2) - 1)))
        print("Parte Real: " + str(eig.real) + " Parte Imaginaria: " + str(eig.imag))
    dmd.plot_eigs(show_axes=True, show_unit_circle=True, filename="images/eigs.png")
    for mode in dmd.modes.T:
        print(mode.real)
        plt.plot(x, mode.real)
        plt.title("Modes")
        plt.savefig("images/modes.png")
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
    '''
main()