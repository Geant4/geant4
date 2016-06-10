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
// $Id: Histo.hh 67994 2013-03-13 11:05:39Z gcosmo $
//
/// \file medical/GammaTherapy/include/Histo.hh
/// \brief Definition of the Histo class
//
#ifndef Histo_h
#define Histo_h 1

//---------------------------------------------------------------------------
//
// ClassName:   Histo
//
// Description: Singleton class to hold Emc geometry parameters.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              Histo::GetPointer() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
// Author:      V.Ivanchenko 27/09/00
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4DataVector.hh"
#include "G4Track.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4RootAnalysisManager;

class Histo
{

public:

  static Histo* GetPointer();

  Histo();

  ~Histo();

  void BeginOfHisto();
  // In this method histogramms are booked

  void EndOfHisto();
  // In this method bookHisto method is called in which histogramms are filled

public: 
  
  void ScoreNewTrack(const G4Track* aTrack);

  void AddPhantomStep(G4double e, G4double r1, G4double z1, 
                      G4double r2, G4double z2,
                      G4double r0, G4double z0);

  void AddPhantomGamma(G4double e, G4double r);

  inline void SetHistoName(const G4String& name)   { fHistName = name; };
  inline void AddStepInTarget()                    { ++fNstepTarget;};

  inline void SetVerbose(G4int val)                { fVerbose = val;};
  inline G4int GetVerbose() const                  { return fVerbose;};

  inline void SetNumberDivZ(G4int val)             { fNBinsZ = val; };
  inline G4int GetNumberDivZ() const               { return fNBinsZ; };

  inline void SetNumberDivR(G4int val)             { fNBinsR = val; };
  inline G4int GetNumberDivR() const               { return fNBinsR; };

  inline void SetNumberDivE(G4int val)             { fNBinsE = val; };

  inline void SetAbsorberZ(G4double val)           { fAbsorberZ = val; };
  inline void SetAbsorberR(G4double val)           { fAbsorberR = val; };
  inline void SetScoreZ(G4double val)              { fScoreZ = val; };

  inline void SetMaxEnergy(G4double val)           { fMaxEnergy = val; };
  inline G4double  GetMaxEnergy() const            { return fMaxEnergy;};

  inline void SetCheckVolume(G4VPhysicalVolume* v) { fCheckVolume = v;};
  inline void SetGasVolume(G4VPhysicalVolume* v)   { fGasVolume = v;};
  inline void SetPhantom(G4VPhysicalVolume* v)     { fPhantom = v; };
  inline void SetTarget1(G4VPhysicalVolume* v)     { fTarget1 = v; };
  inline void SetTarget2(G4VPhysicalVolume* v)     { fTarget2 = v; };
  
private:

  void BookHisto();

  void AddPhantomPhoton(const G4DynamicParticle*);
  void AddTargetPhoton(const G4DynamicParticle*);
  void AddPhantomElectron(const G4DynamicParticle*);
  void AddTargetElectron(const G4DynamicParticle*);

  inline void AddPhoton(const G4DynamicParticle*)  { ++fNgam; };
  inline void AddElectron(const G4DynamicParticle*){ ++fNelec; };
  inline void AddPositron(const G4DynamicParticle*){ ++fNposit; };

  // MEMBERS
  static Histo* fManager;
  G4RootAnalysisManager* fAnalysisManager;

  const G4ParticleDefinition* fGamma;
  const G4ParticleDefinition* fElectron;
  const G4ParticleDefinition* fPositron;

  G4VPhysicalVolume* fCheckVolume;
  G4VPhysicalVolume* fGasVolume;
  G4VPhysicalVolume* fPhantom;
  G4VPhysicalVolume* fTarget1;
  G4VPhysicalVolume* fTarget2;
  G4String fHistName;

  G4int fNHisto;
  G4int fVerbose;
  G4int fNBinsZ;
  G4int fNBinsR;
  G4int fNBinsE;
  G4int fScoreBin;

  G4double fScoreZ;
  G4double fAbsorberZ;
  G4double fAbsorberR;

  G4double fMaxEnergy;

  G4double fStepZ;
  G4double fStepR;
  G4double fStepE;
  G4double fNormZ;
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
  std::vector<G4int> fHisto;
 
};

#endif
