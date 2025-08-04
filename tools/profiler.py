import tomllib
import subprocess
import sys


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


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("error: incorrect usage of profiler.py")
        print("    python profiler.py <simulation-file>.toml")
        exit(1)

    with open(sys.argv[1], "rb") as f:
        sim_prop = tomllib.load(f)

    popen = subprocess.run(
        ["rm", "-rf", "../" + sim_prop["simulation"]["output"]["outputdir"]])
    popen = subprocess.run(
        ["mkdir", "-p", "../" + sim_prop["simulation"]["output"]["outputdir"]])

    runs = get_int_input(":: Enter desired number of profile runs:")

    totaltime = 0.0

    apollo_time_file = "/tmp/apollotime"
    runtime_min = float("inf")
    runtime_max = float("-inf")

    for i in range(0, runs):
        # -P passes a relative path to fix issues with data location and whatnot.
        popen = subprocess.run(
            ["../build/apollo", "-C", "../config/config.toml", "-S", sys.argv[1], "-P", "../"])
        
        runtime = float(open(apollo_time_file).read())
        totaltime += runtime
        runtime_min = min(runtime_min, runtime)
        runtime_max = max(runtime_max, runtime)

    print(":: Total kernel time    : " + str(totaltime))
    print(":: Average kernel time  : " + str(totaltime / runs))
    print(":: Fastest              : " + str(runtime_min))
    print(":: Slowest              : " + str(runtime_max))
