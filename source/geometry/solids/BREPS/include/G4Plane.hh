// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Plane.hh,v 1.1 1999-01-07 16:07:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
