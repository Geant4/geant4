// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FluoTestPrimaryGeneratorAction.hh,v 1.10 2001-10-31 12:33:43 elena Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestPrimaryGeneratorAction_h
#define FluoTestPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class FluoTestDetectorConstruction;
class FluoTestPrimaryGeneratorMessenger;
class FluoTestRunAction;
class FluoTestAnalysisManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

#ifdef G4ANALYSIS_USE    
 
 FluoTestPrimaryGeneratorAction(FluoTestDetectorConstruction*,FluoTestAnalysisManager*);
#else

    FluoTestPrimaryGeneratorAction(FluoTestDetectorConstruction*); 
#endif
   
   ~FluoTestPrimaryGeneratorAction();

  public:
  void GeneratePrimaries(G4Event*);
  void SetRndmFlag(G4String val) { rndmFlag = val;}
  void SetRndmPart (G4String val) { rndmPart = val;}
  void SetRndmVert (G4String val) { beam = val;}
  
private:
  G4ParticleGun*                particleGun;	  //pointer a to G4 service class
  FluoTestDetectorConstruction*    FluoTestDetector;  //pointer to the geometry

  FluoTestPrimaryGeneratorMessenger* gunMessenger; //messenger of this class
  FluoTestRunAction*  runManager;

#ifdef G4ANALYSIS_USE   
    FluoTestAnalysisManager* analysisManager;
#endif


  G4String                      rndmFlag;	  //flag for a random impact point  
  G4String                      rndmPart;//flag for a random particle shot
  G4String                      beam;
  
     
};

#endif


