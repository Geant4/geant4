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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "exrdmEnergyDeposition.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exrdmHisto;

class exrdmAnalysisManager
{

public:
  // With description

  static exrdmAnalysisManager* getInstance();
  static void dispose();

private:

  exrdmAnalysisManager();
  ~exrdmAnalysisManager();

public: // Without description

  void bookHisto();

  void BeginOfRun();
  void EndOfRun(G4int);

  void BeginOfEvent();
  void EndOfEvent();

  void AddParticle(G4double, G4double, G4double, G4double);
  void AddIsotope(G4double, G4double, G4double);
  void AddEnergy(G4double, G4double, G4double);

  void SetVerbose(G4int val) {verbose = val;};
  G4int GetVerbose() const {return verbose;};

  void SetFirstEventToDebug(G4int val) {nEvt1 = val;};
  G4int FirstEventToDebug() const {return nEvt1;};
  void SetLastEventToDebug(G4int val) {nEvt2 = val;};
  G4int LastEventToDebug() const {return nEvt2;};

  void SetMaxEnergyforHisto(G4double val) {histEMax = val;};
  G4double  GetMaxEnergyforHisto() const {return histEMax;};
  void SetMinEnergyforHisto(G4double val) {histEMin = val;};
  G4double  GetMinEnergyforHisto() const {return histEMin;};
  void SetNumBinforHisto(G4int val) {histNBin = val;};
  G4int  GeNumBinforHisto() const {return histNBin;};

  void SetThresholdEnergyforTarget(G4double val) {targetThresE = val;};
  G4double GetThresholdEnergyforTarget () const {return targetThresE;};
  void SetThresholdEnergyforDetector(G4double val) {detectorThresE = val;};
  G4double GetThresholdEnergyforDetector () const {return detectorThresE;};
  void SetPulseWidth(G4double val) {pulseWidth = val;};
  G4double GetPulseWidth () const {return pulseWidth;};

private:

  // MEMBERS
  static exrdmAnalysisManager* fManager;

  G4int verbose;
  G4int nEvt1;
  G4int nEvt2; 

  G4double histEMax;
  G4double histEMin;
  G4int histNBin;
  /*
  G4bool histTarget;
  G4bool histDetector;
  G4bool histCoin;
  G4bool histAntiCD;
  G4bool histAntiCT;
  G4bool histEmission;
  G4bool ntupleEmission;
  G4bool ntupleIsotope;
  */
  G4double targetThresE;
  G4double detectorThresE;
  G4double pulseWidth;

  // energy depositions for an event
  std::vector <exrdmEnergyDeposition> Edepo;
  //
  exrdmHisto*  histo;
  
};

#endif
