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
/// \file radioactivedecay/rdecay02/include/exrdmAnalysisManager.hh
/// \brief Definition of the exrdmAnalysisManager class
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef exrdmAnalysisManager_h
#define exrdmAnalysisManager_h 1

//---------------------------------------------------------------------------
//
// ClassName:   exrdmAnalysisManager
//
// Description: Singleton class to hold analysis parameters and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              exrdmAnalysisManager::GetInstance() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "exrdmEnergyDeposition.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class exrdmHisto;

class exrdmAnalysisManager
{

public:
  // With description

  static exrdmAnalysisManager* GetInstance();
  static void Dispose();

private:

  exrdmAnalysisManager();
  ~exrdmAnalysisManager();

public: // Without description

  void BookHisto();

  void BeginOfRun();
  void EndOfRun(G4int);

  void BeginOfEvent();
  void EndOfEvent();

  void AddParticle(G4double, G4double, G4double, G4double);
  void AddIsotope(G4double, G4double, G4double);
  void AddEnergy(G4double, G4double, G4double);
  void AddDecayProduct(G4double pid,G4int Z, G4int A,
                             G4double energy, G4double time,G4double weight);

  void SetVerbose(G4int val) {fVerbose = val;};
  G4int GetVerbose() const {return fVerbose;};

  void SetFirstEventToDebug(G4int val) {fNEvt1 = val;};
  G4int FirstEventToDebug() const {return fNEvt1;};
  void SetLastEventToDebug(G4int val) {fNEvt2 = val;};
  G4int LastEventToDebug() const {return fNEvt2;};

  void SetMaxEnergyforHisto(G4double val) {fHistEMax = val;};
  G4double  GetMaxEnergyforHisto() const {return fHistEMax;};
  void SetMinEnergyforHisto(G4double val) {fHistEMin = val;};
  G4double  GetMinEnergyforHisto() const {return fHistEMin;};
  void SetNumBinforHisto(G4int val) {fHistNBin = val;};
  G4int  GeNumBinforHisto() const {return fHistNBin;};

  void SetThresholdEnergyforTarget(G4double val) {fTargetThresE = val;};
  G4double GetThresholdEnergyforTarget () const {return fTargetThresE;};
  void SetThresholdEnergyforDetector(G4double val) {fDetectorThresE = val;};
  G4double GetThresholdEnergyforDetector () const {return fDetectorThresE;};
  void SetPulseWidth(G4double val) {fPulseWidth = val;};
  G4double GetPulseWidth () const {return fPulseWidth;};

private:

  // MEMBERS
  static exrdmAnalysisManager* fManager;

  G4int fVerbose;
  G4int fNEvt1;
  G4int fNEvt2; 

  G4double fHistEMax;
  G4double fHistEMin;
  G4int fHistNBin;

  G4double fTargetThresE;
  G4double fDetectorThresE;
  G4double fPulseWidth;

  // energy depositions for an event
  std::vector <exrdmEnergyDeposition> fEdepo;
  //
  exrdmHisto*  fHisto;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
