# Unit Test (proto) Documentation

The core philosophy of this unit test framework is as follows:

- Minimize possibility of an inaccurate test result
- Make the framework as kernel-agnostic as possible
- Handle any possible error (within reason) gracefully

## Logic

The current sample code (in place of a kernel) is trivial to understand. This
allows us to focus on the logic of the unit test. The entry point is in 
`lib/main.c`. The code is simple enough to follow with explanation throughout
as necessary.

## Pedantry

Naturally a level of dogma/pedantry is necessary in this. For the time being
my constraints are:

- C. Data is data. In an effort to minimize the possibility of the compiler
  giving us issues I've written this in C. A sufficiently advanced CPP
  implementation could probably give it a run for it's money, but it would
  take more effort than it's probably worth.
- Kernel agnosticism. The code in lib/ should __NOT__ rely on anything in 
  the kernels.
- Minimize voodoo. There are certainly parts of this that could be made
  "better" with some voodoo. I want this to be maintainable.

## A note on build systems.

At the present moment the GCC compiler is set to use the same options for 
the unit tests as for the apollo source code itself. This is not reflective
of the ultimate aim of these tests. In time we will enforce stricter
IEEE 754 conformance on the testing components to ensure the accuracy of
Apollo's results. However, that is outside the scope of this initial outline.

## Snapshot file data structure.

The first 8 bytes are a validation string (info below).
The next 4 bytes are the signed integer representing the number of params. 
This is the length of both the `void** param` and the `int* size`. The next 
`4 * len` bytes are the size(s). The next `4 * size[i]` bytes (for all `len`
values in `size`. The last 8 bytes are the ending validation string.

## Validation 

A very rudimentary implementation of data validation is implemented in the 
framework. There are a handful of specific circumstances under which it could
fail, as it is intended mainly as a placeholder. The endgame is likely to be 
a combination of the current method, _and_ a hash of the total data (appended
to the end), checked against a master list (to simplify validating many known
datasets).
