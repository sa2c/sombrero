/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include "logger.h"
#include "global.h"

/* 
 * ***********************************************
 * Simple output logging facility
 * ***********************************************
 */

typedef struct _record {
	char *name;
	FILE *file;
	struct _record *next;
} record;

typedef struct _lrecord {
	char *name;
	int level;
	struct _lrecord *next;
} lrecord;

static record *filemap=0; /* list the mappings to files */
static record *files=0; /* list of open files */
static record *default_out=0; /* this is a special list with only one element if stdout is mapped */
static lrecord *levels=0; /* list of verbosity levels of IDs */
static int mapchanged=1; /* to keep track of changes to previous maps */

static int verblevel=0; /* the default verbosity level */


/* find a record with the given file in the list curr */

/* find a record with the given name in the record list *curr */
static record *findname(record *curr, char *name) {

	assert(name!=0);
	while(curr!=0) {
		if(strcmp(name,curr->name)==0)
			return curr;
		curr=curr->next;
	}
	return 0;
}

/* find a record with the given name in the verbosity list */
static lrecord *lfindname(lrecord *curr, char *name) {

	assert(name!=0);
	while(curr!=0) {
		if(strcmp(name,curr->name)==0)
			return curr;
		curr=curr->next;
	}
	return 0;
}

/* create a new record and put it at the beginning of the given list */

/* create a new record and put it at the beginning of the given list */
static void addlrecord(lrecord **list, char* name, int level) {
	lrecord *new;

	assert(name!=0);

	new=malloc(sizeof(*new));
	new->name=malloc(sizeof(char)*(strlen(name)+1));
	strcpy(new->name,name);
	new->level=level;
	new->next=*list;
	*list=new;

	mapchanged=1;

}

/* remove a record from a list */

/* remove a record from a list */


/* reset mappping */

/* link the ID name to the file with name filename */
/* Returns:
 * 0 => success
 * 1 => invalid name
 * 2 => invalid filename
 * 3 => cannot open new file
 */

void logger_setlevel(char *name, int v){
	lrecord *rd;

	if(name==0) { /* set the default verbosity level */
		verblevel=v;
		return;
	}

	rd=lfindname(levels,name);
	if(rd!=0){
		rd->level=v;
		return;
	}

	addlrecord(&levels,name,v);

}

void logger_set_input(input_logger *logger){
  if(logger->def_log_lvl==-1)
    logger->def_log_lvl=10;
  
  logger_setlevel(0,logger->def_log_lvl);


  if(logger->inverter_log_lvl==-1)
    logger->inverter_log_lvl=logger->def_log_lvl;
  else {
    logger_setlevel("INVERTER",logger->inverter_log_lvl);
  }
  if(logger->forcestat_log_lvl==-1)
    logger->forcestat_log_lvl=logger->def_log_lvl;
  else {
    logger_setlevel("FORCE-STAT",logger->forcestat_log_lvl);
  }
}


/* this function reset the verbosity level of name to default */


static void mycpyname(char **dst, char *src){
	*((*dst)++)='[';
	while(*src) {
		*((*dst)++)=*(src++);
	}
	*((*dst)++)=']';
}

static int mycpytonl(char **dst, char **src){
	while(**src) {
		*((*dst)++)=**src;
		if(*((*src)++)=='\n') {
			if(**src=='\0')
				return 0;
			else 
				return 1;
		}
	}
	return 0;
}

static int logger_inactive=0;


void logger_disable() {
  logger_inactive=1;
}

int lprintf(char *name, int level, char *format, ...) {
	va_list args;
	static record *lastrec=0;
	static char lastname[512]={0};
	static FILE *lastfd=0;
	static char buf[1024]; 
	static char alevel[16];
	char *cur=&buf[0];
	int ret;
	lrecord *vrd;
	static int lastvlevel;
	static int newline=1;
	int islast=1;

	if(logger_inactive || name==0) /* no name: print nothing and return and error */
		return -1;

	/* compare current name with last name if map has not changed */
	if(mapchanged || strcmp(name,lastname)!=0) {
		islast=0;
		mapchanged=0;
		lastrec=findname(filemap,name);
		if(lastrec==0){
			lastfd=(default_out==0)?stdout:default_out->file;
		} else {
			lastfd=lastrec->file;
		}
		strcpy(&lastname[0],name);
		vrd=lfindname(levels,name);
		lastvlevel=verblevel;
		if(vrd!=0)
			lastvlevel=vrd->level;
	}

	/* check verbosity level */
	if(lastvlevel<level)
		return 0;

	va_start(args, format);

	sprintf(alevel,"%d",level);
	if(newline) {
		mycpyname(&cur,name);
		mycpyname(&cur,alevel);
	} else if (!islast) {
		*(cur++)='\n';
		mycpyname(&cur,name);
		mycpyname(&cur,alevel);
	}
	while(mycpytonl(&cur,&format)){
		mycpyname(&cur,name);
		mycpyname(&cur,alevel);
	}
	newline=(*(cur-1)=='\n')?1:0;
	*cur='\0';

	ret=vfprintf(lastfd,&buf[0],args);
#ifdef IO_FLUSH
	fflush(lastfd);
#endif

	va_end(args);

	return ret;
	
}



