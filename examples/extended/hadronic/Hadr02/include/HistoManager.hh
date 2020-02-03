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
/// \file hadronic/Hadr02/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
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
// Author:      V.Ivanchenko 03/03/2011
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Histo;
class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VModularPhysicsList;
class G4VPhysicsConstructor;

class HistoManager
{
public:

  static HistoManager* GetPointer();

  ~HistoManager();

  void BookHisto();

  void BeginOfRun();
  void EndOfRun();
  void BeginOfEvent();
  void EndOfEvent();
  void Fill(G4int id, G4double x, G4double w);

  void ScoreNewTrack(const G4Track*);
  void AddTargetStep(const G4Step*);

  void SetVerbose(G4int val);        
  void SetIonPhysics(const G4String&);

  inline void SetTargetMaterial(const G4Material* mat) {fMaterial = mat;};
  inline void SetTargetLength(G4double val)            {fLength  = val;};
  inline void SetNumberOfSlices(G4int val)             {fNSlices = val;};
  inline void SetNumberOfBinsE(G4int val)              {fNBinsE  = val;};

  inline G4double Length()         const               {return fLength;};
  inline G4int    NumberOfSlices() const               {return fNSlices;};

  inline G4int GetVerbose() const                      {return fVerbose;};

  inline void SetDefaultBeamPositionFlag(G4bool f)     {fBeamFlag = f;};
  inline G4bool DefaultBeamPosition() const            {return fBeamFlag;};

  inline void SetMaxEnergyDeposit(G4double val)        {fEdepMax = val;};

  inline void SetPhysicsList(G4VModularPhysicsList* p) {fPhysList = p;};

private:
  HistoManager();

  void Initialise();

  static HistoManager* fManager;

  G4VModularPhysicsList*      fPhysList;
  G4VPhysicsConstructor*      fIonPhysics;
  const G4ParticleDefinition* fPrimaryDef;
  const G4ParticleDefinition* fNeutron;
  const G4Material*           fMaterial;

  G4double fEdepMax;
  G4double fEdepEvt;
  G4double fEdepSum;
  G4double fEdepSum2;
  G4double fLength;
  G4double fAbsZ0;
  G4double fPrimaryKineticEnergy;
 
  G4int fVerbose;
  G4int fNBinsE;
  G4int fNSlices;

  G4int fNevt;
  G4int fNelec;
  G4int fNposit;
  G4int fNgam;
  G4int fNcpions;
  G4int fNpi0;
  G4int fNkaons;
  G4int fNmuons;
  G4int fNions;
  G4int fNdeut;
  G4int fNalpha;
  G4int fNneutron;
  G4int fNproton;
  G4int fNaproton;
  G4int fNstep;
  G4int fNHisto;

  G4bool fBeamFlag;

  Histo* fHisto;
};

#endif
