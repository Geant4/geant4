// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALPrimaryGeneratorAction.hh,v 1.2 2002-10-02 19:40:09 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FCALPrimaryGeneratorAction_h
#define FCALPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class FCALTestbeamSetup;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FCALPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    FCALPrimaryGeneratorAction();    
   ~FCALPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);
  //    void SetRndmFlag(G4String val) { rndmFlag = val;}

  private:
    G4ParticleGun*                particleGun;	  //pointer a to G4 service class
    
  
  private:
  G4int Ievent;
  G4double* X;
  G4double* Y;
  G4double* Z;
  G4double* Cos_X;
  G4double* Cos_Y;
  G4double* Cos_Z;

  G4int Nevent;

};

#endif


