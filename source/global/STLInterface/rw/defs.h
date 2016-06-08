#ifndef __defs
#define __defs

#ifndef G4USE_OLDSTL
#  include <cstdio>
#  include <cstring>
#else
#  include <stdio.h>
#  include <string.h>
#endif

#define RWDEFAULT_CAPACITY (64)

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

const size_t RW_NPOS = ~(size_t)0;
#define CLHEP_MAX_MIN_DEFINED
#include <CLHEP/config/TemplateFunctions.h>

#ifdef G4_HAVE_BOOL
  typedef bool RWBoolean;
#else
  typedef HepBoolean RWBoolean;
#endif

class RWBoundsErr
{
public:

  RWBoundsErr(const char* s,int b=0,int v=0)
    {
      char btmp[80],vtmp[80];
      sprintf(btmp,"%d",b);
      sprintf(vtmp,"%d",v);
      str=new char[strlen(s)+8+strlen(btmp)+7+strlen(vtmp)+1];
      strcpy(str,s);
      strcat(str,": bound:");
      strcat(str,btmp);
      strcat(str," value:");
      strcat(str,vtmp);
    }

  RWBoundsErr(const RWBoundsErr&e)
    {
      str=new char[strlen(e.str)+1];
      strcpy(str,e.str);
    }

  ~RWBoundsErr()
    {
      delete str;
    }

  const char * why() const
    {
      return str;
    }

private:

  char* str;

};
      
class RWGeneralException
{
public:

  RWGeneralException(const char* s)
    {
      str=new char[strlen(s)+1];
      strcpy(str,s);
    }

  RWGeneralException(const RWGeneralException&e)
    {
      str=new char[strlen(e.str)+1];
      strcpy(str,e.str);
    }

  ~RWGeneralException()
    {
      delete str;
    }

  const char * why() const
    {
      return str;
    }

private:

  char* str;

};
  
#define RWTHROW(a) throw a

#endif










