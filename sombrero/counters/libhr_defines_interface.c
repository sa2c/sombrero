#include "libhr_defines_interface.h"

#ifndef TESTCOUNTERS
#include "global.h"    // T,X,Y,Z and borders
#include "suN_types.h" // NF
#endif

float cT() { return T; }
float cX() { return X; }
float cY() { return Y; }
float cZ() { return Z; }
float cNF() { return NF; }
float cTBORDER() { return T_BORDER; }
float cXBORDER() { return X_BORDER; }
float cYBORDER() { return Y_BORDER; }
float cZBORDER() { return Z_BORDER; }
