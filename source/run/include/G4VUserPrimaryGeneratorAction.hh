// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VUserPrimaryGeneratorAction.hh,v 1.1 1999-01-07 16:14:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4VUserPrimaryGeneratorAction_h
#define G4VUserPrimaryGeneratorAction_h 1

class G4Event;

class G4VUserPrimaryGeneratorAction
{
  public:
    G4VUserPrimaryGeneratorAction();
    virtual ~G4VUserPrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event* anEvent) = 0;
};

#endif


