// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PPrimaryVertex.ddl,v 2.0 1998/07/02 16:12:39 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
//


#ifndef G4PPrimaryVertex_h
#define G4PPrimaryVertex_h 1

#include "globals.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PPrimaryVertex 
: public HepPersObj
{
  public:
      G4PPrimaryVertex();
      G4PPrimaryVertex(const G4PrimaryVertex* vertex);
      ~G4PPrimaryVertex();

  private:
      G4double X0;
      G4double Y0;
      G4double Z0;
      G4double T0;
      G4PrimaryParticle * theParticle;

      G4PPrimaryVertex* nextVertex;
      G4int numberOfParticle;
};

#endif

