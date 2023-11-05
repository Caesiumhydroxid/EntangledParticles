# Create a figure and axis
fig, ax = plt.subplots()

# Set the x-axis limits
ax.set_xlim(x_min,x_max)

# Initialize an empty point for the plot
point, = ax.plot([], [])

# Function to initialize the plot
def init():
    point.set_data([], [])
    return point,

# Function to update the plot for each frame
def update(frame):
    x = xs
    y = np.real(sol.y[:,frame][0::2])
    point.set_data(x, y)
    return point,

ani = FuncAnimation(fig, update, frames=len(sol.y[0]), init_func=init, blit=True)
plt.show()
