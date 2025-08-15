# Kernels

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

### Kernel trigger vs preprocessing
A kernel trigger is intended ONLY to interact with module-specific code, 
whereas the preprocessing kernel interacts with other module code
(hydro, mainly).

### Note on Kernel MP
There are two implementations of apollo expected by the build system: apollo, 
and apollo-rocm. apollo-rocm is linked with ROCm, and is passed the 
`__MP_ROCM` definition. This allows an rocm and non-rocm implementation. 
Non-rocm is _required_, rocm is _not_ although if the module is particularly 
computationally expensive, it is desired.

Each kernel module for a set of zones/cells will be passed a single device
to work on.
