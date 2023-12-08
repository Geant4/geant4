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
/// \file hadronic/Hadr01/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
// Description: Singleton class to hold parameters and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              HistoManager::GetPointer() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
// Author:      V.Ivanchenko 27/09/00
//
// Modified:
// 04.06.2006 Adoptation of Hadr01 (V.Ivanchenko)
// 03.10.2006 Add csFlag (V.Ivanchenko)
// 16.11.2006 Add beamFlag (V.Ivanchenko)
//
//----------------------------------------------------------------------------
//

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4Element.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Histo;
class G4Track;
class G4Step;
class G4ParticleDefinition;

class HistoManager
{
public:

  static HistoManager* GetPointer();

private:

  HistoManager();

public: 

  ~HistoManager();

  void BookHisto();

  void BeginOfRun();
  void EndOfRun();
  void BeginOfEvent();
  void EndOfEvent();
  void Fill(G4int id, G4double x, G4double w);

  void ScoreNewTrack(const G4Track*);
  void AddTargetStep(const G4Step*);
  void AddLeakingParticle(const G4Track*);

  void SetVerbose(G4int val);        

  inline void SetTargetRadius(G4double val) { fRadius = val; fR2 = val*val; };
  inline void SetTargetLength(G4double val) { fLength = val; };
  inline void SetNumberOfSlices(G4int val)  { fNSlices = val; };
  inline void SetNumberOfBinsE(G4int val)   { fNBinsE = val; };
  inline void SetDefaultBeamPositionFlag(G4bool f) { fBeamFlag = f; };        
  inline void SetMaxEnergyDeposit(G4double val) { fEdepMax = val; };

  inline G4double Radius() const { return fRadius; };
  inline G4double Length() const { return fLength; };
  inline G4bool   DefaultBeamPosition() const { return fBeamFlag; };
  inline G4int    NumberOfSlices() const { return fNSlices; };
  inline G4int    GetVerbose() const { return fVerbose; };

  inline G4int PrintBertiniXS() const { return fPrintBertiniXS; };
  inline void SetPrintBertiniXS(G4int key) { fPrintBertiniXS = key; };

private:

  static HistoManager* fManager;

  const G4ParticleDefinition* fPrimaryDef;
  const G4ParticleDefinition* fNeutron;

  G4double fR2 = 0.0;
  G4double fRadius;
  G4double fLength;
  G4double fEdepMax;
  G4double fEdepEvt = 0.0;
  G4double fEdepEM = 0.0;
  G4double fEdepPI = 0.0;
  G4double fEdepP = 0.0;
  G4double fEdepSum = 0.0;
  G4double fEdepSum2 = 0.0;
  G4double fAbsZ0 = 0.0;
  G4double fPrimaryKineticEnergy = 0.0;

  G4int fVerbose = 0;
  G4int fNBinsE = 100;
  G4int fNSlices = 300;

  G4int fNevt = 0;
  G4int fNelec = 0;
  G4int fNposit = 0;
  G4int fNgam = 0;
  G4int fNprot_leak = 0;
  G4int fNpiofNleak = 0;
  G4int fNcpions = 0;
  G4int fNpi0 = 0;
  G4int fNkaons = 0;
  G4int fNmuons = 0;
  G4int fNions = 0;
  G4int fNdeut = 0;
  G4int fNalpha = 0;
  G4int fNneutron = 0;
  G4int fNproton = 0;
  G4int fNaproton = 0;
  G4int fNneu_forw = 0;
  G4int fNneu_leak = 0;
  G4int fNneu_back = 0;
  G4int fNstep = 0;
  G4int fNHisto = 28;
  G4int fPrintBertiniXS = -1; // 0 - all

  G4bool fBeamFlag = true;
  G4bool fHistoBooked = false;

  Histo* fHisto;
};

#endif
