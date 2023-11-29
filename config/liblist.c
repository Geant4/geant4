
/*
Given a "libname.map" file on standard input and a list or directory
of .d dependency files liblist produces:
  a) without -l option, a library list ordered according to the libname.map,
     giving a warning if conflicting dependencies are found in the .d files.
  b) with -l option, another libname.map, ordered according firstly to the
     original libname.map file, then reordered according to the dependencies
     found in the .d files.  This option is used for compiling the file
     libname.map from all the dependencies.
  c) with -m <lpath> option, the whole existing libraries list ordered
     according to the libname.map, where libraries are placed in <lpath>.
The .d files are specified in the argument(s).
The libname.map is on standard input.

Usage:
  liblist *.d < libname.map
  liblist -d <ddir> < libname.map
  liblist -l *.d < libname.map
  liblist -ld <ddir> < libname.map
  liblist -l -d <ddir> < libname.map
  liblist -m <lpath> < libname.map
where:
  <ddir>  is a directory name of a directory which is recursively
          searched for dependency files
  <lpath> is the path where libraries are located

Frank Behner, John Allison 13th February 1999.
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>

#define BUFSIZE 1000000
#define TRIGSIZE 1000
#define NLIBMAX 200

extern char *optarg;
extern int optind, opterr, optopt;

char** parsedir(char *directory,int *argc)
{
  DIR *actualdir;
  FILE *actualfile;
  struct dirent *entry; 
  char *buffer=0;
  struct stat status;
  char **targv=0,**ptr,**phelp;
  int len,targc,s;

  /*Open the actual directory*/
  actualdir=opendir(directory);

  if(!actualdir) return targv;

  /*Loop over all entries */
  for(entry=readdir(actualdir);entry!=NULL;entry=readdir(actualdir))
    {
      /* Throw away . and .. */
      if(strcmp(entry->d_name,".")==0 ||
	 strcmp(entry->d_name,"..")==0) continue;
      /* Obtain the status information of that entry */
      if(buffer) free(buffer);
      buffer=(char*) malloc((strlen(directory)+
			     strlen(entry->d_name)+2)*sizeof(char));
      strcpy(buffer,directory);
      strcat(buffer,"/");
      strcat(buffer,entry->d_name);
      s=stat(buffer,&status);
      if(s==0)
	{
	  if(S_ISDIR(status.st_mode))
	    {
	      /* a directory, so we are going recursive*/
	      targc=0;
	      ptr=parsedir(buffer,&targc);
	      if(targc)
		{
		  phelp=targv;
		  targv=(char**) malloc((*argc+targc)*sizeof(char*));
		  memcpy(targv,phelp,*argc*sizeof(char*));
		  memcpy(&targv[*argc],ptr,targc*sizeof(char*));
		  *argc+=targc;
		  free(phelp);
		  free(ptr);
		}
	    }
	  else if(S_ISREG(status.st_mode))
	    {
	      /* a regular file is it .d? */
	      len=strlen(entry->d_name);
	      if(entry->d_name[len-2]=='.' && entry->d_name[len-1]=='d')
		{
		  phelp=targv;
		  targv=(char**) malloc((*argc+1)*sizeof(char*));
		  memcpy(targv,phelp,*argc*sizeof(char*));
		  targv[*argc]=strdup(buffer);
		  (*argc)++;
		  free(phelp);
		}
	    }
	}
      else
	{
	  fprintf
	    (stderr,
	     "  No status - perhaps file %s does not exist.\n",
	     directory);
	  exit(1);
	}
    }

  if(buffer) free(buffer);
  closedir(actualdir);

  return targv;
}
		
 
int main (int argc, char** argv) {

  char static buffer[BUFSIZE],*bufferPtr,workbuf[256];
  char *ptr,*p,**pp,**pp1,**pp2,*directory=0,*libpath=0;
  char **rargv;
  char *libname=0;
  int i,optl=0,optm=0,swapping,c,rargc;
  FILE *fp;

#if defined ( _WIN32 ) || defined ( __CYGWIN__ ) || defined ( __CYGWIN32__ )
  char *ntg4tmp=0,*ntg4tmp1=0;
  int nti; 
#endif

  struct libmap_
    {
      char *lib;      /* Library name, e.g., G4run. */
      char *trigger;  /* Source directory, e.g., source/run/.  */
      int used;       /* True if used by these dependency files.  */
      char **uses;    /* List of library names which this library uses.  */
      struct libmap_ *next;
  };

  struct libmap_ *libmap=0,*libmapPtr=0,*libmapPtr1=0,*libmapPtr2=0,
    *prevPtr1,*prevPtr2,*tmpp,*userLibmapPtr;

  while((c=getopt(argc,argv,"ld: m:"))!=EOF)
    {
      switch(c)
	{
	case 'l':
	  optl=1;
	  break;
	case 'd':
	  directory=strdup(optarg);
	  break;
        case 'm':
          optm=1;
	  libpath=strdup(optarg);
          break;
	}
    }
	  
  /*Adjust parameters after parsing options */

  if(optind<argc)
    {
      rargv=&argv[optind];
      rargc=argc-optind;
    }
  else
    {
      rargv=0;
      rargc=0;
    }

  if(directory)
    {
      if(rargc==0)
	{
	  rargv=parsedir(directory,&rargc);
	}
      else
	{
	  fprintf
	    (stderr,
	     "  ERROR: If you specify a directory don't also specify files\n");
	  exit(1);
	}
    }

  if(optl)fprintf(stderr,"  Reading library name map file...\n");
  while (!feof(stdin))
    {
      /* Get library name... */
      fgets(buffer,BUFSIZE,stdin);
      if(feof(stdin)) break;
      if (strlen(buffer) >= BUFSIZE-1)
      {
         fprintf(stderr,
	   " Internal ERROR: BUFSIZE too small to read library name map file\n");
	 exit(1);
      }
	       /* discarded trailing \n, as gets() was doing */
      if ( buffer[strlen(buffer)-1] == '\n') 
        {   buffer[strlen(buffer)-1]='\0'; }  

      ptr=strtok(buffer,":\n");

      /* Check for duplicate names... */
      for(libmapPtr1=libmap;libmapPtr1;libmapPtr1=libmapPtr1->next)
	{
	  if(strcmp(libmapPtr1->lib,ptr)==0)
	    {
	      fprintf(stderr,"  ERROR: Duplicate library name: %s\n",ptr);
	      fprintf(stderr,
		      "  Perhaps a duplicate subdirectory with"
		      " a GNUmakefile with the same library name.\n"
		      );
	      exit(1);
	    }
	}

      if(libmap)
	{
	  libmapPtr->next=(struct libmap_*) malloc(sizeof(struct libmap_));
	  libmapPtr=libmapPtr->next;
	}
      else /* First time through anchor libmapPtr to libmap. */
	{
	  libmap=(struct libmap_*) malloc(sizeof(struct libmap_));
	  libmapPtr=libmap;
	}
      libmapPtr->next=0;      
      libmapPtr->lib=strdup(ptr);
      libmapPtr->used=0;
      libmapPtr->uses=(char**)calloc(NLIBMAX,sizeof(char*));

      /* If option -l not specified, fill uses list... */
      if(!optl && !optm)
	{
	  pp=libmapPtr->uses;
	  if(ptr)
	    {
	      ptr=strtok(NULL," \n");
	      while (ptr)
		{
		  *pp=strdup(ptr);
		  pp++;
		  ptr=strtok(NULL," \n");
		}
	    }
	}

      if(!optm)
        {
          /* Get directory name... */
          fgets(buffer,BUFSIZE,stdin);
	  if (strlen(buffer) >= BUFSIZE-1)
	  {
             fprintf(stderr,
	       " Internal ERROR: BUFSIZE too small to read directory name\n");
	     exit(1);
	  }
	       /* discarded trailing \n, as gets() was doing */
          if ( buffer[strlen(buffer)-1] == '\n') 
             {   buffer[strlen(buffer)-1]='\0'; }  

          ptr=strtok(buffer,"/");
          if(!ptr)
	    {
	      fprintf(stderr,"  ERROR: \"/\" before \"source\" expected.\n");
	      exit(1);
	    }
          while(ptr&&strcmp (ptr,"source"))ptr=strtok(NULL,"/");
          ptr=strtok(NULL,"/");
          if(!ptr)
	    {
	      fprintf(stderr,"  ERROR: \"source\" expected.\n");
	      exit(1);
	    }
          libmapPtr->trigger=(char*)malloc(TRIGSIZE);
	  if(strlen(ptr)>TRIGSIZE)
          {
            fprintf(stderr,"  ERROR: String overflow for: %s\n", ptr);
            exit(1);
          }
          strcpy(libmapPtr->trigger,ptr);
          ptr=strtok(NULL,"/");
          while(ptr&&strcmp(ptr,"GNUmakefile"))
	    {
	      strcat(libmapPtr->trigger,"/");
	      strcat(libmapPtr->trigger,ptr);
	      ptr=strtok(NULL,"/");
	    }
          if(!ptr)
	    {
	      fprintf
	        (stderr,
	         "  ERROR: \"source/<unique-sub-path>/GNUmakefile\" expected.\n");
	      exit(1);
	    }
	}
    }

  if(optl)fprintf(stderr,"  Reading dependency files...\n");

#if defined ( _WIN32 ) || defined ( __CYGWIN__ ) || defined ( __CYGWIN32__ )
      ntg4tmp=getenv("G4TMP");
      if ( ! ntg4tmp ) 
        {
           fprintf(stderr,"  ERROR: Cannot find environment variable G4TMP\n");
           exit(1);
        }
      ntg4tmp1=strdup(ntg4tmp);
#endif

  for(i=0;i<rargc;i++)
    {
      fp=fopen(rargv[i],"r");
      fgets(buffer,BUFSIZE,fp);

#if defined ( _WIN32 ) || defined ( __CYGWIN__ ) || defined ( __CYGWIN32__ )
      ptr=strchr(ntg4tmp1,':');

      while ( ptr=strchr(buffer,'\\') ) *ptr='/';
 
      while (ntg4tmp1!=NULL &&  (ptr=strstr(buffer,ntg4tmp1))!=NULL )
        {
          for(nti=0;nti<strlen(ntg4tmp1);nti++) ptr[nti]=' ';
        }
#endif
     
      /* Clip target out of dependency file... */
      ptr=strtok(buffer,":");
      
      /* Look for a "user" library... */
      for(libmapPtr=libmap;libmapPtr;libmapPtr=libmapPtr->next)
	{
	  if(strlen(libmapPtr->lib)>256)
          {
            fprintf(stderr,"  ERROR: String overflow for: %s\n", libmapPtr->lib);
            exit(1);
          }
          strcpy(workbuf,libmapPtr->lib);
	  /* Add trailing "/" to distinguish track/ and tracking/, etc. */
	  strcat(workbuf,"/");
	  if(strstr(ptr,workbuf)) break;
	}
      if(libmapPtr)
	{
	  userLibmapPtr=libmapPtr;
	}
      else
	{
	  userLibmapPtr=0;
	}

      if(!optm)
	{
	  /* Look for a "used" library and add it to the "user" uses list... */
	  bufferPtr=strtok(NULL,"\n");  /* Start *after* ":". */
	  if (!bufferPtr) 
	    {
	      fprintf(stderr,"  WARNING: It seems there is nothing after \':\' in dependency file %s.\n", rargv[i]);
	    }
	  else {
	  do
	    {
	      for(libmapPtr=libmap;libmapPtr;libmapPtr=libmapPtr->next)
	        {
		  /* Look for trigger string. */
	          if(strlen(libmapPtr->trigger)>256)
                  {
                    fprintf(stderr,"  ERROR: String overflow for: %s\n", libmapPtr->trigger);
                    exit(1);
                  }
                  strcpy(workbuf,libmapPtr->trigger);
		  strcat(workbuf,"/include");
		  ptr=strstr(bufferPtr,workbuf);
		  if(ptr && (userLibmapPtr != libmapPtr))
		    {
		      libmapPtr->used=1;
		      if(userLibmapPtr)
		        {
			  for(pp=userLibmapPtr->uses;*pp;pp++)
			    {
			      if(strcmp(*pp,libmapPtr->lib)==0)break;
			    }
			  if(!*pp)*pp=libmapPtr->lib;
		        }
		    }
		  /* Also look for library name in case header files are
		     placed in temporary directories under a subdirectory
		     with the same name as the library name.  This can
		     happen with Objectivity which makes header files
		     from .ddl files and places them in a temporary
		     directory. */
	          if(strlen(libmapPtr->lib)>256)
                  {
                    fprintf(stderr,"  ERROR: String overflow for: %s\n", libmapPtr->lib);
                    exit(1);
                  }
		  strcpy(workbuf,libmapPtr->lib);
		  strcat(workbuf,"/");
		  ptr=strstr(bufferPtr,workbuf);
		  if(ptr && (userLibmapPtr != libmapPtr))
		    {
		      libmapPtr->used=1;
		      if(userLibmapPtr)
		        {
			  for(pp=userLibmapPtr->uses;*pp;pp++)
			    {
			      if(strcmp(*pp,libmapPtr->lib)==0)break;
			    }
			  if(!*pp)*pp=libmapPtr->lib;
		        }
		    }
	        }
	      fgets(buffer,BUFSIZE,fp);
	      bufferPtr=buffer;

#if defined ( _WIN32 ) || defined ( __CYGWIN__ ) || defined ( __CYGWIN32__ )
	      while ( ptr=strchr(buffer,'\\') ) *ptr='/';

	      while (ntg4tmp1 &&  (ptr=strstr(buffer,ntg4tmp1)) )
	        {
		  for(nti=0;nti<strlen(ntg4tmp1);nti++) ptr[nti]=' ';
	        }
#endif

	    } while(!feof(fp));
	  fclose(fp);
        }
      }
    }

#if defined ( _WIN32 ) || defined ( __CYGWIN__ ) || defined ( __CYGWIN32__ )
      free(ntg4tmp1);
#endif

  if(optl) /* This option is used for compiling the file libname.map
	      from all the dependencies. */
    {
      fprintf(stderr,"  Checking for circular dependencies...\n");
      for(libmapPtr=libmap;libmapPtr;libmapPtr=libmapPtr->next)
	{
	  for(pp=libmapPtr->uses;*pp;pp++)
	    {
	      for(libmapPtr1=libmap;libmapPtr1!=libmapPtr;
		  libmapPtr1=libmapPtr1->next)
		{
		  if(strcmp(libmapPtr1->lib,*pp)==0)
		    {
		      for(pp1=libmapPtr1->uses;*pp1;pp1++)
			{
			  if(strcmp(*pp1,libmapPtr->lib)==0)break;
			}
		      if(*pp1)
			{
			  fprintf
			    (stderr,
			     "  WARNING: %s and %s use each other.\n",
			     libmapPtr->lib,
			     libmapPtr1->lib);
			}
		    }
		  else 
		    {
		      /* Not right yet...
		      for(pp1=libmapPtr1->uses;*pp1;pp1++)
			{
			  for(libmapPtr0=libmap;libmapPtr0!=libmapPtr1;
			      libmapPtr0=libmapPtr0->next)
			    {
			      if(libmapPtr0==*pp)
				{
				  fprintf
				    (stderr,
				     "  WARNING: triangular dependecy:\n"
				     "  %s uses %s uses %s uses %s.\n",
				     libmapPtr->lib,
				     libmapPtr1->lib,
				     libmapPtr0->lib,
				     libmapPtr->lib);
				}
			    }
			}
		      */
		    }
		}
	    }
	}

      fprintf(stderr,"  Reordering according to dependencies...\n");
      do
	{
	  swapping=0;
	  prevPtr2=0;
	  for(libmapPtr=libmap;libmapPtr;libmapPtr=libmapPtr->next)
	    {
	      for(pp=libmapPtr->uses;*pp;pp++)
		{
		  prevPtr1=0;
		  for(libmapPtr1=libmap;libmapPtr1!=libmapPtr;
		      libmapPtr1=libmapPtr1->next)
		    {
		      if(strcmp(libmapPtr1->lib,*pp)==0)
			{
			  /* Check that 1st doesn't use 2nd... */
			  for(pp1=libmapPtr1->uses;*pp1;pp1++)
			    {
			      if(strcmp(*pp1,libmapPtr->lib)==0)break;
			    }
			  if(!*pp1) /* If not... */
			    {
			      swapping=1;
			      /* Make previous of 1st now point to 2nd... */
			      if(prevPtr1)
				{
				  prevPtr1->next=libmapPtr;
				}
			      else
				{
				  libmap=libmapPtr;
				}
			      /* Make 2nd now point to what 1st did, unless
				 it's adjacent, in which case make it point
				 to 1st itself... */
			      tmpp=libmapPtr->next;
			      if(libmapPtr1->next==libmapPtr)
				{
				  libmapPtr->next=libmapPtr1;
				}
			      else
				{
				  libmapPtr->next=libmapPtr1->next;
				}
			      /* Make previous of 2nd point to 1st, unless
				 it's adjacent, in which case leave it... */
			      if(libmapPtr1->next!=libmapPtr)
				{
				  prevPtr2->next=libmapPtr1;
				}
			      /* Make 1st now point to what 2nd did... */
			      libmapPtr1->next=tmpp;
			      break;
			    }
			}
		      prevPtr1=libmapPtr1;
		    }
		  if(swapping)break;
		}
	      prevPtr2=libmapPtr;
	      if(swapping)break;
	    }
	}while(swapping);

      fprintf(stderr,"  Writing new library map file...\n");
      for(libmapPtr=libmap;libmapPtr;libmapPtr=libmapPtr->next)
	{
	  printf("%s:",libmapPtr->lib);
	  for(pp=libmapPtr->uses;*pp;pp++)
	    {
	      printf(" %s",*pp);
	    }
	  printf("\n");
	  printf("source/%s/GNUmakefile\n",libmapPtr->trigger);
	}
    }
  else if (optm)
    {
      /* create tmp. string libname */ 
      int libname_usable_size=24;
      if ( ! libname ) libname=malloc(libname_usable_size+16);

      /* Write out full library list... */      
      for(libmapPtr=libmap;libmapPtr;libmapPtr=libmapPtr->next)
      {
	if ( strlen(libpath)+strlen(libmapPtr->lib) > libname_usable_size ) {
	  libname_usable_size=(strlen(libpath)+strlen(libmapPtr->lib))*2;
	  free(libname);
	  libname=malloc(libname_usable_size+16);
	}
        /* Check existance of libraries and print out only installed ones */
	  
	
	snprintf(libname, libname_usable_size+16, "%s/lib%s.a", libpath, libmapPtr->lib);
        if (access(libname,R_OK))
	{
	  snprintf(libname, libname_usable_size+16, "%s/lib%s.so", libpath, libmapPtr->lib);
          if (!access(libname,R_OK))
	  {
            printf("-l%s ",libmapPtr->lib);
	  }
          else  /* case MacOS .dylib */
          {
	    snprintf(libname, libname_usable_size+16, "%s/lib%s.dylib", libpath, libmapPtr->lib);
            if (!access(libname,R_OK))
	    {
              printf("-l%s ",libmapPtr->lib);
	    }
          }
	}
        else
	{
          printf("-l%s ",libmapPtr->lib);
	}
	libmapPtr=libmapPtr->next;
      }
    }
  else
    {
      /* Add dependent libraries... */
      for(libmapPtr=libmap;libmapPtr;libmapPtr=libmapPtr->next)
	{
	  if(libmapPtr->used)
	    {
	      for(pp=libmapPtr->uses;*pp;pp++)
		{
		  for(libmapPtr1=libmap;libmapPtr1;libmapPtr1=libmapPtr1->next)
		    {
		      if(strcmp(libmapPtr1->lib,*pp)==0)
			{
			  libmapPtr1->used=1;
			}
		    }
		}	      
	    }
	}

      /* Write out library list... */
      for(libmapPtr=libmap;libmapPtr;libmapPtr=libmapPtr->next)
	{
	  if(libmapPtr->used)
	    {
	      printf("-l%s ",libmapPtr->lib);
	    }
	}
    }
      
  exit(0);
  
}
