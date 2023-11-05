import numpy as np
import matplotlib.pyplot as plt
from project_math import *
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation

x_min = -1e-7
x_max = 1e-7
B0 = 1
xs = np.linspace(x_min,x_max,N)

gaussian = generateGaussianWavePacket(xs,-5e-8,5e-9,9e8)
gaussian2 = generateGaussianWavePacket(xs,0,1e-9,-3e9)
As = generateMagneticVectorPotential(xs,B0,0.5e-8,1e-8)
print(As)
plt.plot(xs,As)
plt.savefig("potential.svg")

gaussianTogether = np.zeros(len(gaussian)*2, dtype=complex)
for i in range(len(gaussian)):
    gaussianTogether[2*i] = gaussian[i]
    #gaussianTogether[2*i+1] = gaussian2[i]

H = hamiltonian(xs,As)

def f(t, phi):
    return 1/(1j*constants.hbar) * (H @ phi)

tmax = 15e-13
sol = solve_ivp(f, [0, tmax], gaussianTogether,  t_eval=np.linspace(0,tmax,100), method='DOP853')


print(sol)
ys = np.abs(sol.y[:,int(0)][0::2])**2

# Create a figure and axis
fig, (ax1,ax2) = plt.subplots(2,1)

# Set the x-axis limits
ax1.set_xlim(x_min,x_max)
ax2.set_xlim(x_min,x_max)

ax1.set_ylim(0,1e8)
ax2.set_ylim(0,1e8)

ax1.title.set_text("Spin Up")
ax2.title.set_text("Spin Down")

# Initialize an empty point for the plot
point1, = ax1.plot(xs, np.real(sol.y[:,int(10)][0::2]))
ax1.plot(xs, As*1e14)
point2, = ax2.plot(xs, np.real(sol.y[:,int(10)][1::2]))
ax2.plot(xs, As*1e14)
#Plot on second y axis


# Function to initialize the plot
def init():
    point1.set_data(xs, ys)
    point2.set_data(xs, ys)
    return point1, point2

# Function to update the plot for each frame
def update(frame):
    x = xs
    print(frame)
    point1.set_data(xs, np.abs(sol.y[:,int(frame)][0::2])**2)
    point2.set_data(xs, np.abs(sol.y[:,int(frame)][1::2])**2)
    return point1, point2

ani = FuncAnimation(fig, update, frames=len(sol.y[0]), init_func=init, blit=True,interval=1/60*1000)
ani.save(f"testSimplN{N}_B{B0}.gif")

plt.show()
