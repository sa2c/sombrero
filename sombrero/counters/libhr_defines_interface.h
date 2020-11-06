#ifndef __LIBHR_DEFINES_INTERFACE_H_
#define __LIBHR_DEFINES_INTERFACE_H_

#ifdef MKPYMOD
def cT() : return T;
def cX() : return X;
def cY() : return Y;
def cZ() : return Z;
def cNF() : return NF;
def cTBORDER() : return T_BORDER;
def cXBORDER() : return X_BORDER;
def cYBORDER() : return Y_BORDER;
def cZBORDER() : return Z_BORDER;
#else
// These cannot be static, they are part of the interface
float cT();
float cX();
float cY();
float cZ();
float cNF();
float cTBORDER();
float cXBORDER();
float cYBORDER();
float cZBORDER();
#endif
#endif // __LIBHR_DEFINES_INTERFACE_H_
