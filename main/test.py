import matplotlib.pyplot as plt
import numpy as np

from pydmd import DMD


def f1(x, t):
    return 1.0 * np.exp(x*1j * t)


def f2(x, t):
    return 2.0 / np.cosh(x) * np.tanh(x) * np.exp(2.8j * t)


x = np.linspace(-5, 5, 128)
t = np.linspace(0, 4 * np.pi, 256)

xgrid, tgrid = np.meshgrid(x, t)

print(len(xgrid))

X1 = f1(xgrid, tgrid)
X2 = f2(xgrid, tgrid)
X = X1
print(X)

titles = ["$f_1(x,t)$", "$f_2(x,t)$", "$f$"]
data = [X1, X2, X]

fig = plt.figure(figsize=(17, 6))
for n, title, d in zip(range(131, 134), titles, data):
    plt.subplot(n)
    plt.pcolor(xgrid, tgrid, d.real)
    plt.title(title)
plt.colorbar()
plt.show()

dmd = DMD(svd_rank=2)
dmd.fit(X.T)

for eig in dmd.eigs:
    print(
        "Eigenvalue {}: distance from unit circle {}".format(
            eig, np.abs(eig.imag ** 2 + eig.real ** 2 - 1)
        )
    )

dmd.plot_eigs(show_axes=True, show_unit_circle=True)

for mode in dmd.modes.T:
    plt.plot(x, mode.real)
    plt.title("Modes")
plt.show()

for dynamic in dmd.dynamics:
    plt.plot(t, dynamic.real)
    plt.title("Dynamics")
plt.show()
