/***************************************************************************\
* Copyright (c)                                                             *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*

int iup_wrk(int site, int dir)
    Up geometry pointer for the active workspace

int idn_wrk(int site, int dir)
    Down geometry pointer for the active workspace

suNg * pu_gauge_wrk(int site, int dir);
    Pointer to the active workspace gauge field poistion and direction element

suNg_field *u_gauge_wrk()
    Pointer to the active workspace gauge field

void reset_wrk_pointers()
    Reset the workspace pointers to point to the default gauge field (u_gauge,iup,idn)

void set_wrk_space(int i)
    Set the workspace pointers to point to the workspace "i" (u_gauge,iup,idn)

void set_wrk_space_and_pointers(int i, suNg_field **g_wrk_out, int **i_up_wrk_out, int **i_dn_wrk_out)
    Set the workspace pointers to point to the workspace "i" (u_gauge,iup,idn) and gives a direct link to the gauge field pointer and to the geometry pointers

int reserve_wrk_space()
    Reserve one workspace and returns the id of the reserved workspace 

int reserve_wrk_space_with_pointers(suNg_field **g_wrk_out, int **i_up_wrk_out, int **i_dn_wrk_out)
    Reserve one workspace and returns the id of the reserved workspace and gives a direct link to the gauge field pointer and to the geometry pointers

void release_wrk_space(int id_release)
    Release the workspace identified by id_release

void free_wrk_space();
    Free all the workspaces
*/

#include "global.h"
#include "suN.h"
#include "memory.h"
#include <stdlib.h>
#include <string.h>
#include "logger.h"

static suNg_field **_g_wrk = NULL;
static int **_iup_wrk;
static int **_idn_wrk;

suNg_field *_g = NULL;
int *_iup;
int *_idn;

static int *_wrk_reserved = NULL;
static int n_alloc = 0;
static int n_reserved = 0;











void free_wrk_space()
{
    int j;
    if (n_alloc != 0)
    {
        for (j = 0; j < n_alloc; j++)
        {
            free_gfield(_g_wrk[j]);
            free(_iup_wrk[j]);
        }

        free(_g_wrk);
        free(_iup_wrk);
        free(_idn_wrk);
        free(_wrk_reserved);
        n_reserved = n_alloc = 0;
    }
}
