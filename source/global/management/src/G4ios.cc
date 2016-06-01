//History 1998 Nov. 3 Masayasu Nagamatu

#include "G4ios.hh"

#ifdef G4STREAM

#include "G4strstreambuf.hh"

G4strstreambuf G4coutbuf;
G4strstreambuf G4cerrbuf;
ostream G4cout(&G4coutbuf);
ostream G4cerr(&G4cerrbuf);

#endif


