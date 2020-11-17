#ifndef _MEMORY_BASE
#define _MEMORY_BASE

#include "fm_defs.h"

#ifndef MKPYMOD
#define REAL_SIZE 8
#endif

#include "libhr_defines_interface.py.h" // provides cT(),cX(),cY(),cZ() and cNF()

_FD(local_sites, cT() * cX() * cY() * cZ());
_FD(local_sites_e, local_sites() / 2);
_FD(halo_sites, 2 * local_sites() *
                    (cTBORDER() / cT() + cXBORDER() / cX() + cYBORDER() / cY() +
                     cZBORDER() / cZ()));
_FD(halo_sites_e, halo_sites() / 2);

_FD(complex_size, 2 * REAL_SIZE);
_FD(vector_size, cNF() * complex_size());
_FD(spinor_size, 4 * vector_size());

#endif
