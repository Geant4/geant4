// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN03PrimaryGeneratorAction.hh,v 1.1 2002-10-01 13:43:56 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef ExN03PrimaryGeneratorAction_h
#define ExN03PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class FCALTestbeamSetup;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ExN03PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ExN03PrimaryGeneratorAction();    
   ~ExN03PrimaryGeneratorAction();

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


