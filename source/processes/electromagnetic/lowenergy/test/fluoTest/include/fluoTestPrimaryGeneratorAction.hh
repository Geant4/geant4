//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef fluoTestPrimaryGeneratorAction_h
#define fluoTestPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class fluoTestDetectorConstruction;
class fluoTestPrimaryGeneratorMessenger;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
   fluoTestPrimaryGeneratorAction(fluoTestDetectorConstruction*);    
   ~fluoTestPrimaryGeneratorAction();

  public:
    
    void GeneratePrimaries(G4Event*);
    void SetRndmFlag(G4String val) { rndmFlag = val;}
    void SetRandomizePrimary (G4String val){ randomizePrimary = val;}   
  private:
    G4ParticleGun*                particleGun;	  //pointer a to G4 service class
    fluoTestDetectorConstruction*    Detector;  //pointer to the geometry
    
    fluoTestPrimaryGeneratorMessenger* gunMessenger; //messenger of this class
    G4String                      rndmFlag;	  //flag for a random impact point       
  G4String                      randomizePrimary; //flag for a random energy of the 
                                                   //particle
   
    G4double momentum;
    G4double sigmaMomentum;
    G4double sigmaAngle;
 
public:
 
    inline void SetMomentum(G4double val) { momentum = val; }
    inline G4double GetMomentum() const { return momentum; }
    void SetSigmaMomentum(G4double);
    void SetSigmaAngle(G4double);
};

#endif


