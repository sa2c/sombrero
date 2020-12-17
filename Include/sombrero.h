#ifndef __SOMBRERO_H_
#define __SOMBRERO_H_

#include "spinor_field.h"

int setup_process_sombrero(int *argc, char ***argv);
int finalize_process_sombrero();

void flat_source(spinor_field *s);

#endif // __SOMBRERO_H_
