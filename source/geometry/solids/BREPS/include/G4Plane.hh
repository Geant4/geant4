// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Plane.hh,v 2.1 1998/10/20 16:31:28 broglia Exp $
// GEANT4 tag $Name: geant4-00 $
//
#ifndef __G4Plane
#define __G4Plane
#include "globals.hh"
class G4Plane 
{
public:
    G4Plane(){a=b=c=d=0;}
    G4double a,b,c,d;
};

#endif
