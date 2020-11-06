#ifndef _MEMORY_BASE
#define _MEMORY_BASE

#include "fm_defs.h"

#ifdef MKPYMOD
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

#include "libhr_defines_interface.h" // provides cT(),cX(),cY(),cZ() and cNF()

_FD( local_sites, cT()*cX()*cY()*cZ());
_FD( local_sites_e, local_sites()/2);
_FD( halo_sites, 2*local_sites()*(1.0/cT()+1.0/cX()+1.0/cY()+1.0/cZ()));
_FD( halo_sites_e, halo_sites()/2);

_FD( complex_size, 2*REAL_SIZE );
_FD( vector_size, cNF()*complex_size() );
_FD( spinor_size, 4*vector_size() );

#endif
