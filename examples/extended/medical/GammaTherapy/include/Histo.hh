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
#include <vector>
#include "G4DynamicParticle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4DataVector.hh"
#include "G4Track.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

namespace AIDA {
 class ITree;
 class ITuple;
 class IHistogram1D;
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
    void AddPositron(const G4DynamicParticle*) {n_posit++;};
    void AddStepInTarget() {n_step_target++;};
    void SetVerbose(G4int val) {verbose = val;};
    G4int GetVerbose() const {return verbose;};
    void SetHistoNumber(G4int val) {nHisto = val;};
    void SetNtuple(G4bool val) {nTuple = val;};

    void SetNumberDivZ(G4int val) {nBinsZ = val; };
    G4int GetNumberDivZ() const    {return nBinsZ;};
    void SetNumberDivR(G4int val) {nBinsR = val; };
    G4int GetNumberDivR() const    {return nBinsR;};
    void SetNumberDivE(G4int val) {nBinsE = val; };

    void SetFirstEventToDebug(G4int val) {nEvt1 = val;};
    G4int FirstEventToDebug() const {return nEvt1;};
    void SetLastEventToDebug(G4int val) {nEvt2 = val;};
    G4int LastEventToDebug() const {return nEvt2;};

    void SetAbsorberZ(G4double val) {absorberZ = val;};
    void SetAbsorberR(G4double val) {absorberR = val;};
    void SetScoreZ(G4double val)    {scoreZ = val;};

    void SetMaxEnergy(G4double val) {maxEnergy = val;};
    G4double  GetMaxEnergy() const {return maxEnergy;};
    void AddEvent() {n_evt++;};
    void AddStep() {n_step++;};
    void SetCheckVolume(G4VPhysicalVolume* v) {checkVolume = v;};
    void SetGasVolume(G4VPhysicalVolume* v) {gasVolume = v;};
    G4VPhysicalVolume* CheckVolume() const {return checkVolume;};
    G4VPhysicalVolume* GasVolume() const {return gasVolume;};
    void SetPhantom(G4VPhysicalVolume* v) {phantom = v;};
    void SetTarget1(G4VPhysicalVolume* v) {target1 = v;};
    void SetTarget2(G4VPhysicalVolume* v) {target2 = v;};

    void AddStep(G4double e, G4double r1, G4double z1, G4double r2, G4double z2,
                             G4double r0, G4double z0);
    void AddGamma(G4double e, G4double r);
    void ScoreNewTrack(const G4Track* aTrack);

private:

  // MEMBERS
    static Histo* fManager;

    G4VPhysicalVolume* checkVolume;
    G4VPhysicalVolume* gasVolume;
    G4VPhysicalVolume* phantom;
    G4VPhysicalVolume* target1;
    G4VPhysicalVolume* target2;
    G4String histName;
    G4String histType;

    std::vector<AIDA::IHistogram1D*> histo;
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
