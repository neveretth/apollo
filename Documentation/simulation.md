# Simulation

Simulation setup for Apollo is done with the simulation.toml file, specified 
by `-S <file>.toml`.

## Options

### simulation

#### simulation.time
- endtime `bool` (the time the simulation runs in seconds)

#### simulation.output
- output `bool` (whether or not output will be generated)
- outputdir `string` (path of the output directory)
- tres `int` (temporal resolution of output, aka: how many frames are printed)

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

#### simulation.neutrino
- use `bool` (is neutrino sim used)

### visualization 
These are only used by the python `runner.py` script.
- save `bool` (should the final fig/animation be saved?)
- saveas `string` (what should it be saved as?)
