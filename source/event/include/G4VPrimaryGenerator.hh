// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPrimaryGenerator.hh,v 2.1 1998/07/12 02:53:57 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef G4VPrimaryGenerator_h
#define G4VPrimaryGenerator_h 1

class G4Event;

class G4VPrimaryGenerator
{
  public:
     G4VPrimaryGenerator();
     virtual ~G4VPrimaryGenerator();

     virtual void GeneratePrimaryVertex(G4Event* evt) = 0;
};

#endif

