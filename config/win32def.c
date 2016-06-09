/* $Id: win32def.c,v 1.3 2004-06-18 16:10:06 gcosmo Exp $ */
/* +---------------------- Copyright notice -------------------------------+ */
/* | Copyright (C) 1995, Guy Barrand, LAL Orsay, (barrand@lal.in2p3.fr)    | */
/* |   Permission to use, copy, modify, and distribute this software       | */
/* |   and its documentation for any purpose and without fee is hereby     | */
/* |   granted, provided that the above copyright notice appear in all     | */
/* |   copies and that both that copyright notice and this permission      | */
/* |   notice appear in supporting documentation.  This software is        | */
/* |   provided "as is" without express or implied warranty.               | */
/* +---------------------- Copyright notice -------------------------------+ */

/*
  This program is a tool to build .def file
 used by lib.exe to build the DLL on Windows.
  It takes the output of a "dumpbin" over an
 archive library and produces the .def file for
 this library. For example :
   DOS> dumpbin /out:tmp /symbols MyLib.arc
   DOS> win32def.exe MyLib < tmp > MyLib.def
 Note that win32def is a standalone program that 
 can be easily reconstructed with :
   DOS> cl /Fowin32def.exe win32def.c
*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

char** GetWords (char*,char*,int*);
/*****************************************************************************/
int main (
 int aArgn 
,char** aArgs
)
/*****************************************************************************/
{
#define MAX_STRING 2048
  char buffer[MAX_STRING+1];
  int length;


  /*EXETYPE WINDOWS\n\ */
  printf("\
LIBRARY %s\n\
EXPORTS\n",aArgs[1]);

  while(1){
    if(fgets(buffer,MAX_STRING,stdin)==NULL) return EXIT_FAILURE;
    /*
      On some system (NT) editors when saving binary files 
      put \r\n at place of \n ; we then look for \r\n.
    */
    length = strlen(buffer);
    if( (length>=2) && (buffer[length-2]=='\r') && (buffer[length-1]=='\n') ) {
      buffer[length-2] = '\0';
      length--;
      length--;
    } else if((length>=1) && (buffer[length-1]=='\n')) {
      buffer[length-1] = '\0';
      length--;
    }
    if(strstr(buffer,"SECT")==NULL) continue;
    if(strstr(buffer,"External")==NULL) continue;
    if(strstr(buffer,"??_")!=NULL) {
      if(strstr(buffer,"operator/=")==NULL) continue;
      /* Keep operator /= */
    }

    {    
      char** words;
      int    wordn;
      words  = GetWords (buffer," ",&wordn);
      if(wordn>=8) {
	int iword = 7;
        int offset = 0;
	if(strcmp(words[4],"()")==0) {
	  iword = 7;
	} else if(strcmp(words[4],"External")==0) {
	  iword = 6;
	}
	if(words[iword][0]=='_') offset = 1;

        if( (iword==6) && (strstr(buffer,"static class")!=NULL) ){
          /* static class objects are not DATA */
	  printf(" %s\n",words[iword]+offset);
	} else if(iword==6) {
          /* DATA */
	  printf(" %s\tDATA\n",words[iword]+offset);
	} else {
          /* code */
	  printf(" %s\n",words[iword]+offset);
	}

      }
      {int count;
      for(count=0;count<wordn;count++) if(words[count]) free(words[count]);
      if(words) free(words);}
    }
    /*printf("%s\n",buffer);*/
  }
  return EXIT_SUCCESS;
}
/*****************************************************************************/
char** GetWords (
 char* a_string 
,char* a_limiter 
,int* a_number
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#define STRDUP(str)  ((str) != NULL ? (strcpy((char*)malloc((unsigned)strlen(str) + 1), str)) : (char*)NULL)
#define STRDEL(str) {if((str)!=NULL) {free(str);str=NULL;}}
  int count;
  char*  string;
  char*  token;
  int    iline;
  char** list  = NULL;
  int    nline = 0;
/*.........................................................................*/
  if(a_number!=NULL) *a_number      = 0;
  if( (a_string==NULL) || (*a_string=='\0') )  return NULL;
  if(a_limiter==NULL) return NULL;

  string = STRDUP(a_string);
  if(string==NULL) return NULL;

  nline = 16;
  list = (char**)malloc(nline*sizeof(char*));
  if(list==NULL) return NULL;
  iline = 0;

  token = string;
  while(1) { 
    char* pos;
    pos = strstr (token,a_limiter);
    if(pos!=NULL) {
      *pos = '\0';
      if(*token!='\0') {
	if(iline>=nline) { 
	  nline    +=16;
	  list      = (char**)realloc(list,nline*sizeof(char*));
	  if(list==NULL) return NULL;
	}
	list[iline]      = token;
	iline++;
      }
      token = pos + strlen(a_limiter);          
    } else { /*last word*/
      if(*token!='\0') {
	if(iline>=nline) {
	  nline    += 16;
	  list      = (char**)realloc(list,nline*sizeof(char*));
	  if(list==NULL) return NULL;
	}
	list[iline]      = token;
	iline++;
      }
      break;
    }
  }
  
  for(count=0;count<iline;count++) list[count] = STRDUP(list[count]);
  STRDEL(string);

  if(iline==0)  {
    if(list) free(list);
    if(a_number!=NULL) *a_number = 0;
    return             NULL;
  } else {
    if(a_number!=NULL) *a_number = iline;
    return             list;
  }
}
