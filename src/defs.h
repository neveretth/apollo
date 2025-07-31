#ifndef __DEFS_H
#define __DEFS_H

#ifndef __REAL_TYPE_DOUBLE
typedef float real_t;
#else // __REAL_TYPE_FLOAT
typedef double real_t;
#endif

// Currently only used by neutrino kernel
#define dt_tol_high 1e-1
#define dt_tol_low 5e-1

#endif
