// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UVHit.hh,v 1.3 2000-02-25 15:59:45 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __G4UV_Hit
#define __G4UV_Hit

#include "globals.hh"

class G4UVHit
{
public:
    G4UVHit(){u=-1;v=-1; next=0;}	
    G4UVHit(G4double u_hit, G4double v_hit){u = u_hit; v = v_hit; next=0;}
    ~G4UVHit(){}

public:
    G4UVHit * next;
    G4double u, v;
};

#endif
