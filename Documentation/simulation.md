# Simulation

Simulation setup for Apollo is done with the simulation.toml file, specified 
by `-S <file>.toml`.

## Options

### simulation

#### simulation.time
- endtime `bool` (the time the simulation runs in seconds)
- tres `int` (temporal resolution, increases clarity of output,
  does not change physical behavior)
- printkerneltime `bool` (print the time to execute the compute kernel(s))

#### simulation.output
- output `bool` (whether or not output will be generated)
- outputdir `string` (path of the output directory)

#### simulation.resolution
- x `int`
- y `int`
- z `int` (this should be set to 1 for 2D simulation)

#### simulation.hydro
- use `bool` (is hydro sim used)
- outputfile `string` (filename for output)
- effect `bool` (is an effect to be used?)

#### simulation.hydroeffect

##### simulation.hydroeffect.temp
- effect `string` (options: radial, gradient, random)
- base `float` (value to initialize the data with)

##### simulation.hydroeffect.density
- effect `string` (options: radial, gradient, random)
- base `float` (value to initialize the data with)

#### simulation.thermo
- use `bool` (is thermo sim used)
- networkfile `string` (location of network file)
- ratefile `string` (location of rate file)

#### simulation.neutrino
- use `bool` (is neutrino sim used)
- opacityfile `string` (location of opacity file)

### visualization 
These are only used by the python `runner.py` script.
- save `bool` (should the final fig/animation be saved?)
- saveas `string` (what should it be saved as?)
