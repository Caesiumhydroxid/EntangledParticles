import numpy as np
import matplotlib.pyplot as plt
from project_math import *
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation

x_min = -1e-7
x_max = 5e-7
B0 = 1
N = 5001
xs = np.linspace(x_min,x_max,N)

gaussian = generateGaussianWavePacket(xs,-5e-8,5e-9,15e8)
As,Bs = generateMagneticVectorPotential(xs,B0,0.5e-8,5e-7)



gaussianTogether = np.zeros(len(gaussian)*2, dtype=complex)
for i in range(len(gaussian)):
    gaussianTogether[2*i] = 1/np.sqrt(2) * gaussian[i]
    gaussianTogether[2*i+1] = 1/np.sqrt(2) * gaussian[i]

delta = xs[1]-xs[0]
H = hamiltonian(xs,delta,As)

def f(t, phi):
    return 1/(1j*constants.hbar) * (H @ phi)

tmax = 50e-13
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
ax1.plot(xs, Bs*1e7)
point2, = ax2.plot(xs, np.real(sol.y[:,int(10)][1::2]))
ax2.plot(xs, Bs*1e7)
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
ani.save(f"bothSpinDirectionsN{N}_B{B0}.gif")

plt.show()
