/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/


#ifndef INPUT_PAR_H
#define INPUT_PAR_H


/* Global or common variables */
typedef struct _input_rlx {
    
    /* random numbers */
    int rlxd_level, rlxd_seed;
    char rlxd_state[256];
    char rlxd_start[256];
    
} input_rlx;

/* Logger global variables */
typedef struct _input_logger {
/* Logger level defined globally */
/* They are defined at imput level */
/* If you need to separate the log level for a channel insert it here */

  int def_log_lvl;
  int inverter_log_lvl;
  int forcestat_log_lvl;
  
} input_logger;


#endif /* INPUT_PAR_H */
