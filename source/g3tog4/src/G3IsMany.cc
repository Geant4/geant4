// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3IsMany.cc,v 1.1 1999-01-07 16:06:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4ios.hh"

#include "globals.hh"
#include "G3toG4.hh"

G4bool G3IsMany(G4String vonly)
{
    if (vonly == "ONLY" )
        {
            return FALSE;
        } else if (vonly == "MANY") {
            return TRUE;
        } else {
            G4cerr << "Unidentified ONLY/MANY tag: " << vonly << endl;
            return FALSE;
        }
}
