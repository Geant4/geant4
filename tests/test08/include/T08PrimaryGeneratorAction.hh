// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08PrimaryGeneratorAction.hh,v 1.1 1999-01-08 16:35:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef T08PrimaryGeneratorAction_h
#define T08PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class T08DetectorConstruction;
class G4ParticleGun;
class G4Event;
class T08PrimaryGeneratorMessenger;

class T08PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    T08PrimaryGeneratorAction(T08DetectorConstruction* myDC);    
    ~T08PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
    void SetRndmFlag(G4String val) { rndmFlag = val;}

  private:
    G4ParticleGun* particleGun;
    T08DetectorConstruction* myDetector;
    
    T08PrimaryGeneratorMessenger* gunMessenger;
    G4String rndmFlag;       
};

#endif


