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
// Rich advanced example for Geant4
// PadHpdPhotoElectricEffect.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef PadHpdPhotoElectricEffect_h
#define PadHpdPhotoElectricEffect_h 1
#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4DynamicParticle.hh"
#include "G4OpticalPhoton.hh" 
#include "G4Electron.hh"
#include "G4Step.hh"
#include "RichTbMaterialParameters.hh"
class RichTbAnalysisManager;
class RichTbRunConfig;

class PadHpdPhotoElectricEffect : public G4VDiscreteProcess {
public:

  PadHpdPhotoElectricEffect(const G4String&,RichTbRunConfig* );
  virtual ~PadHpdPhotoElectricEffect();
  
  G4bool IsApplicable(const G4ParticleDefinition&  );
  // is applicable for optical photon only

  G4double GetMeanFreePath(const G4Track& ,
				 G4double ,
				 G4ForceCondition* condition);
  //returns infinity (max integer possible) . This means the process does
  // not limit the step, but sets the Forced condition for the DoIt to be
  // invoked at every step. But only at the boundary between Hpd quartz 
  // window and the Hpd photocathode any action be taken.

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
				       const G4Step&  aStep);
  //The method for implementing the HpdPhotoelectric effect.


  G4double getHpdQEff(G4int, G4double);
  // To interpolate the QE from QE data.
  G4double getHpdPhElectronKE() 
  {return HpdPhElectronKE; }
  G4double getPhCathodeToSilDetDist()
  {return PhCathodeToSilDetDist; }
  G4double getHpdDemag(G4int hpdnum) {return DemagnificationFactor[hpdnum]; }
 private:
  //  RichTbAnalysisManager* rAnalysisPhy;
  RichTbRunConfig* rConfig;
  G4String PrePhotoElectricVolName;
  G4String PostPhotoElectricVolName;
  G4double HpdPhElectronKE;
  G4double PhCathodeToSilDetDist;
  G4double PSFsigma;
  vector<G4double>DemagnificationFactor;
  vector<G4double>DemagnificationQuadFactor;
  vector<vector<G4double> >HpdQE;
  vector<vector<G4double> >HpdWabin;

};

#endif
