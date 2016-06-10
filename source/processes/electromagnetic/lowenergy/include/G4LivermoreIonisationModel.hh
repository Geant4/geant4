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
// $Id: G4LivermoreIonisationModel.hh 80788 2014-05-12 09:07:49Z gcosmo $
//
// Author: Luciano Pandola
//         on base of G4LowEnergyIonisation developed by A.Forti and V.Ivanchenko
//
// History:
// -----------
// 12 Jan 2009   L. Pandola   1st implementation. Migration from EM process 
//                            to EM model. Physics is unchanged.
// 23 Oct 2009   L. Pandola   remove un-necessary methods to manage atomic 
//                            deexcitation (done by G4VEmModel)
// 01 Jun 2011   V Ivanchenko general cleanup - all old deexcitation code removed
// 04 Jul 2011   L Pandola    removed unused private member
// 
// -------------------------------------------------------------------
//
// Class description:
// Low Energy Electromagnetic Physics, e- ionisation
// with Livermore Model
// -------------------------------------------------------------------

#ifndef G4LIVERMOREIONISATIONMODEL_HH
#define G4LIVERMOREIONISATIONMODEL_HH 1

#include "G4VEmModel.hh"
#include "globals.hh"

class G4eIonisationCrossSectionHandler;
class G4VEnergySpectrum;
class G4ParticleChangeForLoss;
class G4AtomicTransitionManager;

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
		 
  void SetVerboseLevel(G4int vl) {verboseLevel = vl;};
  G4int GetVerboseLevel(){return verboseLevel;};

protected:

  G4ParticleChangeForLoss* fParticleChange;

private:
 
  G4LivermoreIonisationModel & operator=(const G4LivermoreIonisationModel &right);
  G4LivermoreIonisationModel(const G4LivermoreIonisationModel&);

  //Intrinsic energy limits of the model: cannot be extended by the parent process
  G4double fIntrinsicLowEnergyLimit;
  G4double fIntrinsicHighEnergyLimit;

  G4bool isInitialised;
 
  G4int verboseLevel;
 
  G4eIonisationCrossSectionHandler* crossSectionHandler;
  G4VEnergySpectrum* energySpectrum;

  G4AtomicTransitionManager* transitionManager;
};

#endif

