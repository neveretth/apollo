# Kernels

NOTE: this concept is not fully implemented, and this is an evolving spec.

## What is a kernel?
A kernel contains primary computation code that is executed at some stage of 
Apollo's runtime. They are contained in the kernel directory (in src/). 

## Example
Let's look at the hydrodynamic kernel as an example. The computation code
exists in src/kernel/hydro, within which are a handful of files. The kernel.c
file holds the code, and kernel.h holds the corresponding prototypes that 
will need to be accessed by the main integration loop. kernel.hip holds HIP 
kernels (thought this is not expected to be implemented at the moment), and
finally the Makefile provides instructions for building the kernels.

## Accessing from the main integration loop.
_The following is not implemented yet..._

The main integration loop will activate a trigger corresponding to each kernel
set activated (hydro, thermonuclear, etc) in the simulation.toml file. This 
trigger will call the computational kernels outlined in the module.

NOTE: This is just the outline, expect this to change before release!
