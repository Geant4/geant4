// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02PrimaryGeneratorAction.hh,v 1.3 2000-12-04 16:24:05 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
 
#ifndef ExN02PrimaryGeneratorAction_h
#define ExN02PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class ExN02DetectorConstruction;
class G4ParticleGun;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
 
class ExN02PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ExN02PrimaryGeneratorAction(ExN02DetectorConstruction*);    
   ~ExN02PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun* particleGun;
    ExN02DetectorConstruction* myDetector;
};

#endif


