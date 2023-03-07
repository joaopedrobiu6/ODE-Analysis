import numpy as np
from pydmd import DMD
import scipy as sp
from scipy.integrate import solve_ivp, odeint
import matplotlib.pyplot as plt

delta = 0.7
epsilon = 0.2


def f(t, r):
    omega = r[0]
    theta = r[1]
    return np.array([- (delta + epsilon*np.cos(t))* np.sin(theta), omega])

time = np.linspace(0, 100, 1000)
init_r = [0.1, 0]

results = solve_ivp(f, (0, 100), init_r, method='RK45', t_eval=time, rtol=1e-8)

plt.plot(results.t, results.y[1])
plt.plot(results.t, results.y[0])

x = results.y[0]

tgrid, xgrid = np.meshgrid(results.t, results.y[0])

dmd = DMD(svd_rank=2)
dmd.fit(x.T)

for eig in dmd.eigs:
    print('Eigenvalue {}: distance from unit circle {}'.format(eig, np.abs(eig.imag**2+eig.real**2 - 1)))

dmd.plot_eigs(show_axes=True, show_unit_circle=True)


plt.show()

