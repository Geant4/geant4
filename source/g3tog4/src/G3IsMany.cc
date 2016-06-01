// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3IsMany.cc,v 2.1 1998/07/13 16:50:26 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
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
