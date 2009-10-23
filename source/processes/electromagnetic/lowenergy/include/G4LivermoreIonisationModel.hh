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
// $Id: G4LivermoreIonisationModel.hh,v 1.3 2009-10-23 09:28:37 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 12 Jan 2009   L. Pandola   1st implementation. Migration from EM process 
//                            to EM model. Physics is unchanged.
// 23 Oct 2009   L. Pandola   remove un-necessary methods to manage atomic 
//                            deexcitation (done by G4VEmModel)
// 
// -------------------------------------------------------------------
//
// Class description:
// Low Energy Electromagnetic Physics, e- ionisation
// with Livermore Model
// -------------------------------------------------------------------

#ifndef G4LIVERMOREIONISATIONMODEL_HH
#define G4LIVERMOREIONISATIONMODEL_HH 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4eIonisationCrossSectionHandler.hh"
#include "G4VEnergySpectrum.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicDeexcitation.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4Material;
class G4ShellVacancy;

class G4LivermoreIonisationModel : public G4VEmModel 
{

public:
  
  G4LivermoreIonisationModel(const G4ParticleDefinition* p=0,
			 const G4String& processName = "LowEnergyIoni");
  
  virtual ~G4LivermoreIonisationModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  
  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                              G4double kinEnergy,
                                              G4double Z,
                                              G4double A=0,
                                              G4double cut=0,
                                              G4double emax=DBL_MAX);
					 
  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);
				   
  virtual G4double ComputeDEDXPerVolume(const G4Material*,
                               const G4ParticleDefinition*,
                               G4double kineticEnergy,
                               G4double cutEnergy);
		

  virtual void SampleDeexcitationAlongStep(const G4Material*,
                                           const G4Track&,
                                           G4double& eloss);

  // min cut in kinetic energy allowed by the model
  virtual G4double MinEnergyCut(const G4ParticleDefinition*,
                                const G4MaterialCutsCouple*);
		 
  void SetVerboseLevel(G4int vl) {verboseLevel = vl;};
  G4int GetVerboseLevel(){return verboseLevel;};

  void ActivateAuger(G4bool);

protected:
  G4ParticleChangeForLoss* fParticleChange;

private:
 
  G4LivermoreIonisationModel & operator=(const G4LivermoreIonisationModel &right);
  G4LivermoreIonisationModel(const G4LivermoreIonisationModel&);

  void InitialiseFluorescence();

  //Intrinsic energy limits of the model: cannot be extended by the parent process
  G4double fIntrinsicLowEnergyLimit;
  G4double fIntrinsicHighEnergyLimit;
  G4int fNBinEnergyLoss;

  G4bool isInitialised;
 
  G4int verboseLevel;
 
  G4eIonisationCrossSectionHandler* crossSectionHandler;
  G4VEnergySpectrum* energySpectrum;
  G4ShellVacancy* shellVacancy;

  G4AtomicDeexcitation deexcitationManager;
  const G4AtomicTransitionManager* transitionManager;

};

#endif

