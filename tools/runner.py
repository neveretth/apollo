import tomllib
import subprocess
import matplotlib.pyplot as plt
import matplotlib.animation as aplt
import pandas as pd
import numpy as np
import sys

fig, ax = plt.subplots()


# Display animated figure to screen (generated from list of figs "ims")
def display_fig(ims, sim_prop):
    ani = aplt.ArtistAnimation(
        fig, ims, interval=30, blit=True, repeat_delay=0)
    if (sim_prop["visualization"]["save"]):
        ani.save("export/" + sim_prop["visualization"]["saveas"], dpi=400)
    plt.show()


# Generate graph from linear data.
def graph_linear(sim_prop):
    datafile = "../" + sim_prop["simulation"]["output"]["outputdir"] + "/" + \
        sim_prop["simulation"]["output"]["outputfile"]

    data = pd.read_csv(datafile, delimiter=' ', header=None)
    ims = []
    grid = data.iloc[0].to_numpy()
    im, = ax.plot(grid)
    ims.append([im])
    tres = sim_prop["simulation"]["time"]["tres"]

    ax.set_xlabel("Zone")
    ax.set_ylabel("Temperature")

    for i in range(1, tres):
        grid = data.iloc[i].to_numpy()
        im, = ax.plot(grid, color="Black")
        ims.append([im])

    display_fig(ims, sim_prop)


# Generate graph from 2D data.
def graph_flat(sim_prop):
    size = int(sim_prop["simulation"]["resolution"]["x"])

    datafile = "../" + sim_prop["simulation"]["output"]["outputdir"] + "/" + \
        sim_prop["simulation"]["output"]["outputfile"]

    data = pd.read_csv(datafile, delimiter=' ', header=None)
    ims = []
    grid = data.iloc[0].to_numpy()
    grid = np.reshape(grid, (-1, size))
    im = ax.imshow(grid, cmap='hot',
                   interpolation='none', animated=True)
    clim = im.properties()['clim']
    ims.append([im])

    ax.set_title("Color = Temp")
    ax.set_xlabel("Zone x")
    ax.set_ylabel("Zone y")

    for i in range(1, sim_prop["simulation"]["time"]["tres"]):
        grid = data.iloc[i].to_numpy()
        grid = np.reshape(grid, (-1, size))
        im = ax.imshow(grid, cmap='hot',
                       interpolation='none', animated=True, clim=clim)
        ims.append([im])

    display_fig(ims, sim_prop)


# Generate graph from 3D data.
def graph_rt(sim_prop):
    sizex = int(sim_prop["simulation"]["resolution"]["x"])
    sizey = int(sim_prop["simulation"]["resolution"]["y"])
    sizez = int(sim_prop["simulation"]["resolution"]["z"])

    datafile = "../" + sim_prop["simulation"]["output"]["outputdir"] + "/" + \
        sim_prop["simulation"]["output"]["outputfile"]

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

    for i in range(1, sim_prop["simulation"]["time"]["tres"]):
        rtgrid = np.zeros((sizex, sizey, sizez))
        grid = data.iloc[i].to_numpy()
        grid = np.reshape(grid, (-1, sizex*sizey))
        for j in range(sizex):
            rtgrid[j] = np.reshape(grid[j], (-1, sizey))

        ax = plt.axes(projection="3d")
        ax.scatter(x, y, z, c=rtgrid, marker=".",
                   cmap=plt.hot(), vmin=8e9, vmax=9e9)
        im = ax
        ims.append([im])

    display_fig(ims, sim_prop)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("error: incorrect usage of runner.py")
        print("    python runner.py <simulation-file>.toml")
        exit(1)

    with open(sys.argv[1], "rb") as f:
        sim_prop = tomllib.load(f)

    popen = subprocess.run(
        ["rm", "-rf", "../" + sim_prop["simulation"]["output"]["outputdir"]])
    popen = subprocess.run(
        ["mkdir", "-p", "../" + sim_prop["simulation"]["output"]["outputdir"]])

    # -P passes a relative path to fix issues with data location and whatnot.
    popen = subprocess.run(
        ["../build/apollo", "-C", "../config/config.toml", "-S", sys.argv[1], "-P", "../"])

    if sim_prop["simulation"]["resolution"]["y"] == 1:
        graph_linear(sim_prop)
    elif sim_prop["simulation"]["resolution"]["z"] == 1:
        graph_flat(sim_prop)
    else:
        graph_rt(sim_prop)
