import matplotlib.pyplot as plt
import matplotlib.animation as aplt
import pandas as pd
import numpy as np
import sys

fig, ax = plt.subplots()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: python heatmap.py datafile <squaregriddim>")
        exit(1)

    size = int(sys.argv[2])

    data = pd.read_csv(sys.argv[1], delimiter=' ', header=None)
    ims = []

    x_ = np.linspace(0, 1, size)
    y_ = np.linspace(0, 1, size)
    z_ = np.linspace(0, 1, size)

    x, y, z = np.meshgrid(x_, y_, z_)

    rtgrid = np.zeros((size, size, size))
    grid = data.iloc[0].to_numpy()
    grid = np.reshape(grid, (-1, size*size))
    for j in range(size):
        rtgrid[j] = np.reshape(grid[j], (-1, size))

    ax = plt.axes(projection="3d")
    ax.scatter(x, y, z, c=rtgrid, marker=".",
               cmap=plt.hot(), vmin=8e9, vmax=14e9)
    im = ax
    ims.append([im])

    for i in range(1, 1001):
        rtgrid = np.zeros((size, size, size))
        grid = data.iloc[i].to_numpy()
        grid = np.reshape(grid, (-1, size*size))
        for j in range(size):
            rtgrid[j] = np.reshape(grid[j], (-1, size))

        ax = plt.axes(projection="3d")
        ax.scatter(x, y, z, c=rtgrid, marker=".",
                   cmap=plt.hot(), vmin=8e9, vmax=14e9)
        im = ax
        ims.append([im])

    ani = aplt.ArtistAnimation(
        fig, ims, interval=20, blit=False, repeat_delay=100)
    ani.save("out.mp4", dpi=400)
