import numpy as np
import matplotlib.pyplot as plt

k = 1.0
m = 1./3.

A = np.matrix([[0,    1],
               [-k/m, 0]])
def F(t):
    return np.matrix([[0],
                      [3 * np.cos(t)]])

nt = 0
max_nt = 10000
dt = 0.001
X = np.zeros([2, max_nt])
print A * X[:,1].reshape(-1, 1)
while (nt < max_nt - 1):
    t = nt * dt
    foo = dt * (A * X[:, nt].reshape(-1, 1) + F(t)) + X[:, nt].reshape(-1, 1)
    X[0, nt + 1] = foo[0]
    X[1, nt + 1] = foo[1]
    nt = nt + 1
t = np.linspace(0, dt * max_nt, X.shape[1])

# Load libMesh solution
data = np.loadtxt("../build/solution.txt", comments='#')
print data.shape

plt.figure()
plt.plot(t, X[0, :],
         data[:, 0], data[:, 1])
plt.legend(['Python', 'libMesh'])
# plt.xlim([0, 2 * np.pi])
plt.ylim([-3, 3])

# plt.figure()
# plt.plot(data[:, 0], data[:, 2],
#          data[:, 0], np.sin(data[:, 0]), 'ro')
plt.show()

