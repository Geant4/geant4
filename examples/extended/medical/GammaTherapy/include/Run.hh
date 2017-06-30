//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm5/include/Run.hh
/// \brief Definition of the Run class
//
// $Id: Run.hh 71375 2013-06-14 07:39:33Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4DynamicParticle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4DataVector.hh"
#include "g4root.hh"

#include "globals.hh"

class DetectorConstruction;
class PrimaryGeneratorAction;
class HistoManager;
class G4ParticleDefinition;
class G4Track;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
 public:
  Run(DetectorConstruction*, PrimaryGeneratorAction*, HistoManager*);
  ~Run();

 public:
    
   virtual void Merge(const G4Run*);
   void EndOfRun();   
  
  void ScoreNewTrack(const G4Track* aTrack);
  
  void AddPhantomStep(G4double e, G4double r1, G4double z1, 
                      G4double r2, G4double z2,
                      G4double r0, G4double z0);
  
  void AddPhantomGamma(G4double e, G4double r);
  
  bool GetVerbose() const                          { return fVerbose; }
  inline void AddStepInTarget()                    { ++fNstepTarget;};
  
private:
  
  DetectorConstruction*  fDetector;
  HistoManager* fHistoMgr;
  G4AnalysisManager* fAnalysisManager;

  std::vector<G4int> fHistoId;
  G4int fNHisto;

  void AddPhantomPhoton(const G4DynamicParticle*);
  void AddTargetPhoton(const G4DynamicParticle*);
  void AddPhantomElectron(const G4DynamicParticle*);
  void AddTargetElectron(const G4DynamicParticle*);

  inline void AddPhoton()  { ++fNgam; };
  inline void AddElectron(){ ++fNelec; };
  inline void AddPositron(){ ++fNposit; };

  // Parameters retrived from DetectorConstructor
  const G4ParticleDefinition* fGamma;
  const G4ParticleDefinition* fElectron;
  const G4ParticleDefinition* fPositron;

  const G4VPhysicalVolume* fCheckVolume;
  const G4VPhysicalVolume* fGasVolume;
  const G4VPhysicalVolume* fPhantom;
  const G4VPhysicalVolume* fTarget1;
  const G4VPhysicalVolume* fTarget2;

  G4int fNBinsR;
  G4int fNBinsZ;
  G4int fNBinsE;
  G4int fScoreBin;

  G4double fScoreZ;
  G4double fAbsorberZ;
  G4double fAbsorberR;
  G4double fMaxEnergy;

  G4double fStepZ;
  G4double fStepR;
  G4double fStepE;
  //  G4double fNormZ;

  // Local histogramming parameters
  G4bool fVerbose;
  G4double fSumR;

  G4int fNevt;
  G4int fNelec;
  G4int fNposit;
  G4int fNgam;
  G4int fNstep;
  G4int fNgamPh;
  G4int fNgamTar;
  G4int fNeTar;
  G4int fNePh;
  G4int fNstepTarget;

  G4DataVector fVolumeR;
  G4DataVector fGammaE;
  G4DataVector fEdep;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

