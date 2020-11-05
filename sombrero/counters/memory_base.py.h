#ifndef _MEMORY_BASE
#define _MEMORY_BASE

#include "fc_defs.h"

#ifdef PYTHON
// just to have these declared.
REAL_SIZE = 8
NF = 4
T = 8
X = 8
Y = 8
Z = 8
#else
#define REAL_SIZE 8
#endif

_FD( local_sites, T*X*Y*Z);
_FD( local_sites_e, local_sites()/2);
_FD( halo_sites, 2*local_sites()*(1.0/T+1.0/X+1.0/Y+1.0/Z));
_FD( halo_sites_e, halo_sites()/2);

_FD( complex_size, 2*REAL_SIZE );
_FD( vector_size, NF*complex_size() );
_FD( spinor_size, 4*vector_size() );

#endif
