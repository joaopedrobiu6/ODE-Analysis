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


leitura = read("data1.txt")
print(leitura)
