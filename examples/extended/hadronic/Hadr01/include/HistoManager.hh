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
// $Id: HistoManager.hh,v 1.1 2006-06-02 19:00:00 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef HistoManager_h
#define HistoManager_h 1

//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
// Description: Singleton class to hold Emc parameters and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              HistoManager::GetPointer() static method.
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
#include "G4ThreeVector.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Histo;

class HistoManager
{

public:
  // With description

  static HistoManager* GetPointer();

private:

  HistoManager();

public: // Without description

  ~HistoManager();

  void bookHisto();

  void BeginOfRun();
  void EndOfRun();

  void BeginOfEvent();
  void EndOfEvent();

  void ScoreNewTrack(const G4Track* aTrack);
  void AddEnergy(G4double edep, G4double step, const G4ThreeVector& v);

  void AddElectron(const G4DynamicParticle*);
  void AddPhoton(const G4DynamicParticle*);
  void AddPositron(const G4DynamicParticle*);
  void AddNeutron(const G4DynamicParticle*);
  void AddProton(const G4DynamicParticle*);
  void AddHadron(const G4DynamicParticle*);
  void AddIon(const G4DynamicParticle*, const G4ThreeVector& pos);
  void AddLeakingNeutron(G4double e, const G4ThreeVector& pos);
  void AddLeakingHadron(G4double e, const G4ThreeVector& pos, 
			const G4ParticleDefinition* p);
  void SetEndPoint(const G4ThreeVector& v);
  void SetAbsRadius(G4double val) {absRadius = val;};
  void SetAbsWidth(G4double val)  {absWidth = val;};
  void SetGapWidth(G4double val)  {gapWidth = val;};
  void SetNumberOfAbs(G4int val)  {n_abs = val;};
  void SetNumBinsE(G4int val)     {nBinsE = val;};
  void SetNtuple(G4bool val)      {nTuple = val;};

  G4double AbsRadius()   const {return absRadius;};
  G4double AbsWidth()    const {return absWidth;};
  G4double GapWidth()    const {return gapWidth;};
  G4int    NumberOfAbs() const {return n_abs;};

  G4double GetTrackLength() const   {return trackLength;};
  void ResetTrackLength()           {trackLength = 0.0;};
  void AddTrackLength(G4double x)   {trackLength += x;};
  void SetVerbose(G4int val);
  G4int GetVerbose() const          {return verbose;};

  void SetFirstEventToDebug(G4int val) {nEvt1 = val;};
  G4int FirstEventToDebug() const {return nEvt1;};
  void SetLastEventToDebug(G4int val) {nEvt2 = val;};
  G4int LastEventToDebug() const {return nEvt2;};

  void SetMaxNeutronEnergy(G4double val) {maxNeutronEnergy = val; maxHadronEnergy=val;};
  void SetMaxElectronEnergy(G4double val) {maxElectronEnergy = val;};
  void SetMaxGammaEnergy(G4double val) {maxGammaEnergy = val;};
  void SetNumBinsEnergy(G4int val) {nBinsE = val;};
  void SetNumBinsXY(G4int val) {n_XY = val;};
  void SetNumBinBragg(G4int val) {nBinXY = val;};

  void AddStep(const G4ParticleDefinition* pd, G4double e) {
    n_step++; currentDef = pd; currentKinEnergy = e;};
  G4double CurrentKinEnergy() {return currentKinEnergy;};
  const G4ParticleDefinition* CurrentDefinition() {return currentDef;};

private:

  void IsotopeStudy();

  // MEMBERS
  static HistoManager* fManager;
  const G4ParticleDefinition* primaryDef;
  const G4ParticleDefinition* currentDef;
  G4double currentKinEnergy;
 
  G4int verbose;
  G4int nEvt1;
  G4int nEvt2;
  G4int nBinsE;
  G4int nBinXY;
  G4int n_XY;

  G4double beamEnergy;
  G4double maxEnergy;
  G4double trackLength;
  G4double absRadius;
  G4double absWidth;
  G4double absLength;
  G4double gapWidth;
  G4double absZ0;
  G4double endX2;
  G4double endY2;
  G4double endZ;
  G4double endZ2;
  G4double Z0;
  G4double A0;
  G4double maxNeutronEnergy;
  G4double maxHadronEnergy;
  G4double maxElectronEnergy;
  G4double maxGammaEnergy;
  G4double primaryKineticEnergy;

  G4int n_abs;
  G4int n_evt;
  G4int n_elec;
  G4int n_posit;
  G4int n_gam;
  G4int n_prot_leak;
  G4int n_aprot_leak;
  G4int n_pos_hadr;
  G4int n_neg_hadr;
  G4int n_ions;
  G4int n_neutron;
  G4int n_neu_forw;
  G4int n_neu_leak;
  G4int n_neu_back;
  G4int n_step;
  G4int n_ion;
  G4int nHisto;
  G4int nHisto1;

  G4bool nTuple;

  Histo* histo;
};

#endif
