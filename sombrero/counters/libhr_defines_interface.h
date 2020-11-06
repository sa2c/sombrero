#ifndef __LIBHR_DEFINES_INTERFACE_H_
#define __LIBHR_DEFINES_INTERFACE_H_

#ifdef MKPYMOD
def cT() : return T;
def cX() : return X;
def cY() : return Y;
def cZ() : return Z;
def cNF() : return NF;
#else
// These cannot be static, they are part of the interface
float cT();
float cX();
float cY();
float cZ();
float cNF();
#endif
#endif // __LIBHR_DEFINES_INTERFACE_H_
