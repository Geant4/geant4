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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

namespace AIDA {
  class ITree;
  class ITuple;
  class IHistogram1D;
  class IAnalysisFactory;
}

class Histo
{

public:
  // With description

  static Histo* GetPointer();

  Histo();

  ~Histo();

  void BeginOfHisto();
  // In this method histogramms are booked

  void EndOfHisto();
  // In this method bookHisto method is called in which histogramms are filled

public: // Without description
  
  void SetHistoName(const G4String& name) {histName = name;};
  void SetHistoType(const G4String& type) {histType = type;};
  void bookHisto();
  void SaveToTuple(const G4String&, G4double);
  void SaveToTuple(const G4String&, G4double, G4double);
  void SaveEvent();

  G4double GetTrackLength() const {return trackLength;};
  void ResetTrackLength() {trackLength = 0.0, trackAbs = true;};
  void SetTrackOutAbsorber() {trackAbs = false;};
  G4bool GetTrackInAbsorber() const {return trackAbs;};
  void AddTrackLength(G4double x)   {trackLength += x;};

  void AddDeltaElectron(const G4DynamicParticle*);
  void AddPhoton(const G4DynamicParticle*);
  void AddPhantomPhoton(const G4DynamicParticle*);
  void AddTargetPhoton(const G4DynamicParticle*);
  void AddPhantomElectron(const G4DynamicParticle*);
  void AddTargetElectron(const G4DynamicParticle*);
  inline void AddPositron(const G4DynamicParticle*) { ++n_posit;};
  inline void AddStepInTarget() { ++n_step_target;};

  inline void SetVerbose(G4int val) {verbose = val;};
  inline G4int GetVerbose() const {return verbose;};

  inline void SetHistoNumber(G4int val) {nHisto = val;};
  inline void SetNtuple(G4bool val) {nTuple = val;};

  inline void SetNumberDivZ(G4int val) {nBinsZ = val; };
  inline G4int GetNumberDivZ() const    {return nBinsZ;};
  inline void SetNumberDivR(G4int val) {nBinsR = val; };
  inline G4int GetNumberDivR() const    {return nBinsR;};
  inline void SetNumberDivE(G4int val) {nBinsE = val; };

  inline void SetFirstEventToDebug(G4int val) {nEvt1 = val;};
  inline G4int FirstEventToDebug() const {return nEvt1;};
  inline void SetLastEventToDebug(G4int val) {nEvt2 = val;};
  inline G4int LastEventToDebug() const {return nEvt2;};

  inline void SetAbsorberZ(G4double val) {absorberZ = val;};
  inline void SetAbsorberR(G4double val) {absorberR = val;};
  inline void SetScoreZ(G4double val)    {scoreZ = val;};

  inline void SetMaxEnergy(G4double val) {maxEnergy = val;};
  inline G4double  GetMaxEnergy() const {return maxEnergy;};

  inline void AddEvent() { ++n_evt; };
  inline void AddStep()  { ++n_step; };

  inline void SetCheckVolume(G4VPhysicalVolume* v) {checkVolume = v;};
  inline void SetGasVolume(G4VPhysicalVolume* v) {gasVolume = v;};
  inline G4VPhysicalVolume* CheckVolume() const {return checkVolume;};
  inline G4VPhysicalVolume* GasVolume() const {return gasVolume;};

  inline void SetPhantom(G4VPhysicalVolume* v) {phantom = v;};
  inline void SetTarget1(G4VPhysicalVolume* v) {target1 = v;};
  inline void SetTarget2(G4VPhysicalVolume* v) {target2 = v;};
  
  void AddStep(G4double e, G4double r1, G4double z1, G4double r2, G4double z2,
                             G4double r0, G4double z0);
  void AddGamma(G4double e, G4double r);
  void ScoreNewTrack(const G4Track* aTrack);

private:

  // MEMBERS
  static Histo* fManager;

  const G4ParticleDefinition* gamma;
  const G4ParticleDefinition* electron;
  const G4ParticleDefinition* positron;
  const G4ParticleDefinition* neutron;

  G4VPhysicalVolume* checkVolume;
  G4VPhysicalVolume* gasVolume;
  G4VPhysicalVolume* phantom;
  G4VPhysicalVolume* target1;
  G4VPhysicalVolume* target2;
  G4String histName;
  G4String histType;

  std::vector<AIDA::IHistogram1D*> histo;
  AIDA::IAnalysisFactory* af;  
  AIDA::ITuple* ntup;
  AIDA::ITree* tree;
  G4int nHisto;
  G4int nHisto1;
  G4int verbose;
  G4int nBinsZ;
  G4int nBinsR;
  G4int nBinsE;
  G4int nScoreBin;
  G4int nEvt1;
  G4int nEvt2;

  G4double absorberZ;
  G4double stepZ;
  G4double scoreZ;
  G4double absorberR;
  G4double stepR;
  G4double maxEnergy;
  G4double stepE;
  G4double normZ;
  G4double sumR;

  G4double trackLength;
  G4bool trackAbs;        // Track is in absorber
  G4int n_evt;
  G4int n_elec;
  G4int n_posit;
  G4int n_gam;
  G4int n_step;
  G4int n_gam_ph;
  G4int n_gam_tar;
  G4int n_e_tar;
  G4int n_e_ph;
  G4int n_step_target;
  G4int n_neutron;
  G4bool nTuple;
  G4DataVector   volumeR;
  G4DataVector   gammaE;

};

#endif
