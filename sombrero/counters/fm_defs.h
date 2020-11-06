#ifndef _FM_DEFS
#define _FM_DEFS

#ifdef MKPYMOD
// for one-line function definition
#define _FD(name, body) def name() : return (body)
#else
// for one-line function definition
#define _FD(name, body)                                                        \
  static float name() { return (body); }
#endif

#endif
