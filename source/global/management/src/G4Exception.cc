// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Exception.cc,v 1.8 1999-11-29 10:17:46 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// G4Exception
//
// Global error function prints string to G4cerr, and aborts
// program
//
// History:
// 30.06.95 P.Kent

#include "G4ios.hh"
#include <stdlib.h>
#include "g4rw/cstring.h"

void G4Exception(const char* s)
{
    if (s)
	{
	    G4cerr << s << G4endl;
	}

   G4cerr << G4endl << "*** G4Exception: Aborting execution ***" << G4endl;
   abort();
}

void G4Exception(G4std::string s)
{
  G4Exception(s.c_str());
}

void G4Exception(G4String s)
{
  G4Exception(s.c_str());
}
