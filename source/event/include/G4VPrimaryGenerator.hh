// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPrimaryGenerator.hh,v 1.1 1999-01-07 16:06:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

