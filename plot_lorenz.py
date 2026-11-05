import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("lorenz_points.csv")

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

# Clean up plot
NUM_TRAJECTORIES = 10
POINTS_PER_TRAJECTORY = len(data) // NUM_TRAJECTORIES

for i in range(NUM_TRAJECTORIES):
    start = i * POINTS_PER_TRAJECTORY
    end = (i + 1) * POINTS_PER_TRAJECTORY

    traj = data.iloc[start:end:20]  # Downsample each trajectory

    ax.plot(traj["x"], traj["y"], traj["z"], linewidth=0.5)

ax.set_title("Lorenz Attractor")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.show()
