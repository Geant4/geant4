// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VUserPrimaryGeneratorAction.hh,v 1.2 1999-11-01 03:12:03 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4VUserPrimaryGeneratorAction_h
#define G4VUserPrimaryGeneratorAction_h 1

class G4Event;

// class description:
//
//  This is the abstract base class of the user's mandatory action class
// for primary vertex/particle generation. This class has only one pure
// virtual method GeneratePrimaries() which is invoked from G4RunManager
// during the event loop.
//  Note that this class is NOT intended for generating primary vertex/particle
// by itself. This class should 
//  - have one or more G4VPrimaryGenerator concrete classes such as G4ParticleGun 
//  - set/change properties of generator(s)
//  - pass G4Event object so that the generator(s) can generate primaries.
//

class G4VUserPrimaryGeneratorAction
{
  public:
    G4VUserPrimaryGeneratorAction();
    virtual ~G4VUserPrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event* anEvent) = 0;
};

#endif


