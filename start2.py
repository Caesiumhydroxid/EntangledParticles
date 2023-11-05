import numpy as np
import matplotlib.pyplot as plt
from project_math import *
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation


x_min = -1e-7
x_max = 1e-7
B0 = 1
xs = np.linspace(x_min,x_max,N)
delta = xs[1]-xs[0]

gaussianA = generateGaussianWavePacket(xs,-5e-8,5e-9,9e8)
gaussianB = generateGaussianWavePacket(xs,0,5e-9,-9e8)
As = generateMagneticVectorPotential(xs,B0,5e-10)

gaussianAUp = np.zeros(N*N*2, dtype=complex)
gaussianADown = np.zeros(N*N*2, dtype=complex)
gaussianBUp = np.zeros(N*N*2, dtype=complex)
gaussianBDown = np.zeros(N*N*2, dtype=complex)

for z1 in range(N):
    for z2 in range(N):
        gaussianAUp  [mapZ12StoLin(z1,z2,0,N)] = 1/np.sqrt(N)*gaussianA[z1]
        gaussianADown[mapZ12StoLin(z1,z2,1,N)] = 1/np.sqrt(N)*gaussianA[z1]
        gaussianBUp  [mapZ12StoLin(z1,z2,0,N)] = 1/np.sqrt(N)*gaussianB[z2]
        gaussianBDown[mapZ12StoLin(z1,z2,1,N)] = 1/np.sqrt(N)*gaussianB[z2]


H = hamiltonian2Particles(xs,delta, As)


def f(t, phi):
    return 1/(1j*constants.hbar) * (H @ phi)

tmax = 15e-13
timesteps = 100

startCond = gaussianAUp

sol = solve_ivp(f, [0, tmax], startCond,  t_eval=np.linspace(0,tmax,timesteps), method='DOP853')
print(sol)

solutionUpA   = np.zeros((timesteps,N),dtype=complex)
solutionDownA = np.zeros((timesteps,N),dtype=complex)
solutionUpB   = np.zeros((timesteps,N),dtype=complex)
solutionDownB = np.zeros((timesteps,N),dtype=complex)

for i in range(timesteps):
    for z in range(N):
        solutionUpA  [i][z] = np.sum(np.power(np.abs([(sol.y[mapZ12StoLin(z,x,0,N)][i]) for x in range(N)]),2))
        solutionDownA[i][z] = np.sum(np.power(np.abs([(sol.y[mapZ12StoLin(z,x,1,N)][i]) for x in range(N)]),2))

        solutionUpB  [i][z] = np.sum(np.power(np.abs([(sol.y[mapZ12StoLin(x,z,0,N)][i]) for x in range(N)]),2))
        solutionDownB[i][z] = np.sum(np.power(np.abs([(sol.y[mapZ12StoLin(x,z,1,N)][i]) for x in range(N)]),2))


ys = solutionUpA[0]
print(np.trapz(ys,xs))

# Create a figure and axis
fig, ((ax1,ax3),(ax2,ax4)) = plt.subplots(2,2)

# Set the x-axis limits
ax1.set_xlim(x_min,x_max)
ax2.set_xlim(x_min,x_max)
ax3.set_xlim(x_min,x_max)
ax4.set_xlim(x_min,x_max)

ax1.set_ylim(0,1e8)
ax2.set_ylim(0,1e8)
ax3.set_ylim(0,1e8)
ax4.set_ylim(0,1e8)

# Initialize an empty point for the plot
point1, = ax1.plot(xs, solutionUpA[0])
point2, = ax2.plot(xs, solutionDownA[0])
point3, = ax3.plot(xs, solutionUpB[0])
point4, = ax4.plot(xs, solutionDownB[0])
ax1.plot(xs, As*1e14)
ax2.plot(xs, As*1e14)
ax3.plot(xs, As*1e14)
ax4.plot(xs, As*1e14)

# Function to initialize the plot
def init():
    point1.set_data(xs, ys)
    point2.set_data(xs, ys)
    point3.set_data(xs, ys)
    point4.set_data(xs, ys)
    return point1, point2, point3, point4

# Function to update the plot for each frame
def update(frame):
    x = xs
    point1.set_data(xs, solutionUpA[frame])
    point2.set_data(xs, (solutionDownA[frame]))
    point3.set_data(xs, solutionUpB[frame])
    point4.set_data(xs, (solutionDownB[frame]))
    return point1, point2, point3, point4


ani = FuncAnimation(fig, update, frames=len(sol.y[0]), init_func=init, blit=True,interval=1/60*1000)
ani.save(f"testTwoParticle{N}_B{B0}.gif")
plt.show()