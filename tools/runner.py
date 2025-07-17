import os
import tomllib
import subprocess
import matplotlib.pyplot as plt
import matplotlib.animation as aplt
import pandas as pd
import numpy as np
import sys

fig, ax = plt.subplots()


def graph_flat(sim_prop):  # 2D graph
    size = int(sim_prop["simulation"]["resolution"]["x"])

    datafile = sim_prop["simulation"]["output"]["outputdir"] + "/" + \
        sim_prop["simulation"]["hydro"]["outputfile"]

    print(datafile)

    data = pd.read_csv(datafile, delimiter=' ', header=None)
    ims = []
    grid = data.iloc[0].to_numpy()
    grid = np.reshape(grid, (-1, size))
    im = ax.imshow(grid, cmap='hot',
                   interpolation='none', animated=True)
    clim = im.properties()['clim']
    ims.append([im])

    for i in range(1, sim_prop["simulation"]["output"]["tres"]):
        grid = data.iloc[i].to_numpy()
        grid = np.reshape(grid, (-1, size))
        im = ax.imshow(grid, cmap='hot',
                       interpolation='none', animated=True, clim=clim)
        ims.append([im])

    ani = aplt.ArtistAnimation(
        fig, ims, interval=25, blit=True, repeat_delay=10)
    if (sim_prop["visualization"]["save"]):
        ani.save("export/" + sim_prop["visualization"]["saveas"], dpi=400)
    plt.show()


def graph_rt(sim_prop):  # 3D graph
    sizex = int(sim_prop["simulation"]["resolution"]["x"])
    sizey = int(sim_prop["simulation"]["resolution"]["y"])
    sizez = int(sim_prop["simulation"]["resolution"]["z"])

    datafile = sim_prop["simulation"]["output"]["outputdir"] + "/" + \
        sim_prop["simulation"]["hydro"]["outputfile"]

    data = pd.read_csv(datafile, delimiter=' ', header=None)
    ims = []

    x_ = np.linspace(0, 1, sizex)
    y_ = np.linspace(0, 1, sizey)
    z_ = np.linspace(0, 1, sizez)

    x, y, z = np.meshgrid(x_, y_, z_)

    rtgrid = np.zeros((sizex, sizey, sizez))
    grid = data.iloc[0].to_numpy()
    grid = np.reshape(grid, (-1, sizex*sizey))
    for j in range(sizex):
        rtgrid[j] = np.reshape(grid[j], (-1, sizey))

    ax = plt.axes(projection="3d")
    ax.scatter(x, y, z, c=rtgrid, marker=".",
               cmap=plt.hot(), vmin=8e9, vmax=14e9)
    im = ax
    ims.append([im])

    for i in range(1, sim_prop["simulation"]["output"]["tres"]):
        rtgrid = np.zeros((sizex, sizey, sizez))
        grid = data.iloc[i].to_numpy()
        grid = np.reshape(grid, (-1, sizex*sizey))
        for j in range(sizex):
            rtgrid[j] = np.reshape(grid[j], (-1, sizey))

        ax = plt.axes(projection="3d")
        ax.scatter(x, y, z, c=rtgrid, marker=".",
                   cmap=plt.hot(), vmin=8e9, vmax=14e9)
        im = ax
        ims.append([im])

    ani = aplt.ArtistAnimation(
        fig, ims, interval=20, blit=False, repeat_delay=100)
    if (sim_prop["visualization"]["save"]):
        ani.save("export/" + sim_prop["visualization"]["saveas"], dpi=400)
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("error: incorrect usage of runner.py")
        print("    python runner.py <simulation-file>.toml")
        exit(1)

    with open(sys.argv[1], "rb") as f:
        sim_prop = tomllib.load(f)

    if not os.path.exists(sim_prop["simulation"]["output"]["outputdir"]):
        os.mkdir(sim_prop["simulation"]["output"]["outputdir"])

    popen = subprocess.run(
        ["../build/apollo", "-C", "../config/config.toml", "-S", sys.argv[1]])

    if sim_prop["simulation"]["resolution"]["z"] == 1:
        graph_flat(sim_prop)
    else:
        graph_rt(sim_prop)
