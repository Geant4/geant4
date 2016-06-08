#ifndef __rstream
#define __rstream

#include "rw/defs.h"

#ifndef G4USE_OLDSTL
#include <iostream>
#else
#include <iostream.h>
#endif

inline istream& rwEatwhite(istream& stream) {return stream >> ws;}

#endif
