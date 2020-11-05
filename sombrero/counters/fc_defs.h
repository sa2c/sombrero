#ifndef _FC_DEFS
#define _FC_DEFS

#ifdef PYTHON
// for one-line function definition
#define _FD(name, body) def name(): return (body)
#else
// for one-line function definition
#define _FD(name, body) static float name(){ return (body);}
#endif

#endif
