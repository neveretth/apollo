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
    grid = data.iloc[0].to_numpy()
    grid = np.reshape(grid, (-1, size))
    im = ax.imshow(grid, cmap='hot',
                   interpolation='none', animated=True)
    clim = im.properties()['clim']
    ims.append([im])

    for i in range(1, 101):
        grid = data.iloc[i].to_numpy()
        grid = np.reshape(grid, (-1, size))
        im = ax.imshow(grid, cmap='hot',
                       interpolation='none', animated=True, clim=clim)
        ims.append([im])

    ani = aplt.ArtistAnimation(
        fig, ims, interval=25, blit=True, repeat_delay=10)
    plt.show()
    ani.save("out.mp4")
