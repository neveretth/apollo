# Simulation

Simulation setup for Apollo is done with the simulation.toml file, specified 
by `-S <file>.toml`.

## TOML
TOML (Tom's Obvious Markup Language) is a simple user-facing config language
designed to map to a hash table.
```TOML
[foo]
bar = OxDEADBEEF

[foobar]
[foobar.foo]
bar = "APOLLO"

```
This would be represented as a string by:
```
foo.bar
# and
foobar.foo.bar
```
The last value ("bar" in the example) after is the only part of the
expression, the rest is the section header.
```
[foo] <-- Section header
bar = "VAL" <-- expression

foobar.foo.bar
---------- ___
     |       L this is the "left = right"
     L this is the [section header]
```

Finally, the documentation uses the following:
- The header is the section header
- The values in the table are the left hand of the expression


## Simulation

### simulation.time

| Value | Type | Description |
|-------|------|-------------|
| __endtime__ | `float` | The total runtime of the simulation. |
| __tres__ | `int` | The temporal resolution (data output). |
| __printkerneltime__ | `bool` | Print intermediate kernel time. |
| __initdt__ | `float` | The initial timestep of the hydrodynamic code. |

### simulation.output

| Value | Type | Description |
|-------|------|-------------|
| __output__ | `bool` | Write output to files. |
| __outputdir__ | `string` | Output directory path. |
| __temp__ | `bool` | Output temp. |
| __density__ | `bool` | Output density. |
| __entropy__ | `bool` | Output entropy. |

### simulation.resolution

| Value | Type | Description |
|-------|------|-------------|
| __x__ | `int` | X dimension size. |
| __y__ | `int` | Y dimension size. (1 for linear) |
| __z__ | `int` | Z dimension size. (1 for 2D) |

### simulation.hydro

| Value | Type | Description |
|-------|------|-------------|
| __use__ | `bool` | Trigger hydro compute kernel. |
| __effect__ | `bool` | Use hydro effect. |
| __ttc__ | `float` | Thermal transfer coefficient. |
| __volume__ | `float` | Volume of star. (not implemented) |

### simulation.hydroeffect
#### simulation.hydroeffect.temp

| Value | Type | Description |
|-------|------|-------------|
| __base__ | `float` | Base temperature value |
| __effect__ | `string` | random, gradient, radial1, radial2, radial3 |

#### simulation.hydroeffect.density

| Value | Type | Description |
|-------|------|-------------|
| __base__ | `float` | Base density value |
| __effect__ | `string` | random, gradient, radial1, radial2, radial3 |

### simulation.thermo

| Value | Type | Description |
|-------|------|-------------|
| __use__ | `bool` | Trigger thermo compute kernel. |
| __networkfile__ | `string` | Filepath of network data. |
| __ratefile__ | `string` | Filepath of rate data. |

### simulation.neutrino

| Value | Type | Description |
|-------|------|-------------|
| __use__ | `bool` | Trigger neutrino compute kernel. |
| __opacityfile__ | `string` | Filepath of neutrino opacity data. |
