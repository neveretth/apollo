# Build System

## Building APOLLO
For the non-rocm version:
```
make
# OR
make apollo
```
for the (very incomplete) rocm-enabled version:
```
make apollo-rocm
```

## Makefile structure
Make is the build system for Apollo. The basic structure is: each directory
in which make needs to perform actions contains a Makefile that performs
said actions.

Say you have three directories: the root directory (root), (root)/src, and
(root)/src/module. There would exist a Makefile in each directory, handling 
compilation for the files within each directory. Each file includes the
neccessary Makefile(s) from only the next level directory: so (root)/Makefile
would include (root)/src/Makefile, and (root)/src/Makefile would include 
(root)/src/module/Makefile.
