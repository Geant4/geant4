#ifndef _Definitions_H
#define _Definitions_H

#define MF_EXIT_FAILURE 1  
#define MF_EXIT_SUCCESS 0 

#ifndef NEW
#define NEW new
#endif

#ifdef NO_EXCEPTION_HANDLING
  #include <stdlib.h>
  class Error;
  extern void userErrorMessage(const Error& e);
  #define Throw(t) userErrorMessage(t); exit(MF_EXIT_FAILURE)
#else
  #define Throw(t) throw t
#endif

#endif // _Definitions_H

