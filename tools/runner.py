import tomllib
import subprocess
import matplotlib.pyplot as plt
import matplotlib.animation as aplt
import pandas as pd
import numpy as np
import sys
# Argument parsing
import argparse

parser = argparse.ArgumentParser(
                prog='ApolloVis',
                description='What does ApolloVis do?',
                epilog='')
parser.add_argument('--simulation-file', nargs=1, help='simulation run file', required=True)
parser.add_argument('--run-simulation', help='run new Apollo simulation', action='store_true')




#fig, ax = plt.subplots()

def display_menu():
    menu = """
    ===========================
         Simulation Options
    ===========================
    1. Temperature
    2. Density
    3. Entropy
    0. Quit
    ===========================
    """
    print(menu)
    

def get_bool_input(msg):
    while 1:
        val = input(msg + " ")
        try:
            val = str(val)
            if val[0] == 'y':
                return True
            elif val[0] == 'n':
                return False
            print(":: invalid input? try again.")
        except:
            print(":: invalid input? try again.")

    return False


def get_int_input(msg):
    while 1:
        val = input(msg + " ")
        try:
            val = int(val)
            return val
        except:
            print(":: invalid input? try again.")

    return 0


# Display animated figure to screen (generated from list of figs "ims")
def display_fig(ims, sim_prop):
    ani = aplt.ArtistAnimation(
        fig, ims, interval=30, blit=True, repeat_delay=0)
    plt.show(block=False)
    if get_bool_input(":: save figure? (y/n)"):
        saveas = input(":: save as: ")
        print(":: saving fig to export/" + saveas + ".mp4")
        ani.save("export/" + saveas + ".mp4", dpi=400)


# Generate graph from linear data.
def graph_linear(sim_prop, datafile, label):

    data = pd.read_csv(datafile, delimiter=' ', header=None)
    ims = []
    grid = data.iloc[0].to_numpy()
    im, = ax.plot(grid, color="Black")
    ims.append([im])
    tres = sim_prop["simulation"]["time"]["tres"]

    ax.set_xlabel("Zone")
    ax.set_ylabel(label)

    for i in range(1, tres+1):
        grid = data.iloc[i].to_numpy()
        im, = ax.plot(grid, color="Black")
        ims.append([im])

    display_fig(ims, sim_prop)


# Generate graph from 2D data.
def graph_flat(sim_prop, datafile, label):
    size = int(sim_prop["simulation"]["resolution"]["x"])

    data = pd.read_csv(datafile, delimiter=' ', header=None)

    # Shitty hack to find clim
    clim_min = float("inf")
    clim_max = float("-inf")
    for i in range(0, sim_prop["simulation"]["time"]["tres"]):
        grid = data.iloc[i].to_numpy()
        grid = np.reshape(grid, (-1, size))
        clim_min = min(np.min(grid), clim_min)
        clim_max = max(np.max(grid), clim_max)

    ims = []
    grid = data.iloc[0].to_numpy()
    grid = np.reshape(grid, (-1, size))
    im = ax.imshow(grid, cmap='hot',
                   interpolation='none', animated=True)
    im.set_clim(clim_min, clim_max)
    ims.append([im])

    ax.set_title("Color = " + label)
    ax.set_xlabel("Zone x")
    ax.set_ylabel("Zone y")

    for i in range(1, sim_prop["simulation"]["time"]["tres"] + 1):
        grid = data.iloc[i].to_numpy()
        grid = np.reshape(grid, (-1, size))
        im = ax.imshow(grid, cmap='hot',
                       interpolation='none', animated=True)
        im.set_clim(clim_min, clim_max)
        ims.append([im])

    display_fig(ims, sim_prop)


# Generate graph from 3D data.
def graph_rt(sim_prop, datafile, label):
    print("here")
    sizex = int(sim_prop["simulation"]["resolution"]["x"])
    sizey = int(sim_prop["simulation"]["resolution"]["y"])
    sizez = int(sim_prop["simulation"]["resolution"]["z"])

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

    for i in range(1, sim_prop["simulation"]["time"]["tres"] + 1):
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


def gen_graph(sim_prop, datafile, label):
    if sim_prop["simulation"]["resolution"]["y"] == 1:
        graph_linear(sim_prop, datafile, label)
    elif sim_prop["simulation"]["resolution"]["z"] == 1:
        graph_flat(sim_prop, datafile, label)
    else:
        graph_rt(sim_prop, datafile, label)


if __name__ == "__main__":
    
    args = parser.parse_args()
    print(args)
    simulation_file = args.simulation_file[0]    


    with open(simulation_file, "rb") as f:
        sim_prop = tomllib.load(f)
    
    
    datafile_temp = "../" + \
        sim_prop["simulation"]["output"]["outputdir"] + "/temp.out"
    datafile_density = "../" + \
        sim_prop["simulation"]["output"]["outputdir"] + "/density.out"
    datafile_entropy = "../" + \
        sim_prop["simulation"]["output"]["outputdir"] + "/entropy.out"

    option_map = {
        1: ("Temperature", datafile_temp),
        2: ("Density", datafile_density),
        3: ("Entropy", datafile_entropy),
        }
    
    if args.run_simulation:
        popen = subprocess.run(
            ["rm", "-rf", "../" + sim_prop["simulation"]["output"]["outputdir"]])
        popen = subprocess.run(
            ["mkdir", "-p", "../" + sim_prop["simulation"]["output"]["outputdir"]])
        # -P passes a relative path to fix issues with data location and whatnot.
        popen = subprocess.run(
            ["../build/apollo", "-C", "../config/config.toml", "-S", simulation_file, "-P", "../"])

    

    while True:
        display_menu()
        choice = int(input("Select an option [0–3]: ") or "0")
        if choice == 0:
            print("Exiting the program.")
            break
        
        if choice in option_map:
            label, datafile = option_map[choice]
            fig, ax = plt.subplots()
            gen_graph(sim_prop, datafile, label)
        else:
            print(f"Invalid input: {choice}. Please enter a number between 0 and 3.")
            
            


    
    
