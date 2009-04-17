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
// $Id: G4PenelopeBremsstrahlungModel.hh,v 1.2 2009-04-17 10:29:20 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 05 Dec 2008   L. Pandola   1st implementation. Migration from EM process 
//                            to EM model. Physics is unchanged.
//
// -------------------------------------------------------------------
//
// Class description:
// Low Energy Electromagnetic Physics, e+ and e- bremsstrahlung
// with Penelope Model
// -------------------------------------------------------------------

#ifndef G4PENELOPEBREMSSTRAHLUNGMODEL_HH
#define G4PENELOPEBREMSSTRAHLUNGMODEL_HH 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4VCrossSectionHandler.hh"
#include "G4PhysicsLogVector.hh"
#include "G4AtomicDeexcitation.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4Material;
class G4VEnergySpectrum;
class G4PenelopeBremsstrahlungAngular;
class G4PenelopeBremsstrahlungContinuous;

class G4PenelopeBremsstrahlungModel : public G4VEmModel 
{

public:
  
  G4PenelopeBremsstrahlungModel(const G4ParticleDefinition* p=0,
			 const G4String& processName ="PenelopeBrem");
  
  virtual ~G4PenelopeBremsstrahlungModel();

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
				 
  // min cut in kinetic energy allowed by the model
  virtual G4double MinEnergyCut(const G4ParticleDefinition*,
                                const G4MaterialCutsCouple*);

  void SetVerbosityLevel(G4int lev){verboseLevel = lev;};
  G4int GetVerbosityLevel(){return verboseLevel;};


protected:
  G4ParticleChangeForLoss* fParticleChange;

private:
 
  G4PenelopeBremsstrahlungModel & operator=(const G4PenelopeBremsstrahlungModel &right);
  G4PenelopeBremsstrahlungModel(const G4PenelopeBremsstrahlungModel&);

 
  //Intrinsic energy limits of the model: cannot be extended by the parent 
  // process
  G4double fIntrinsicLowEnergyLimit;
  G4double fIntrinsicHighEnergyLimit;

  G4int verboseLevel;

  G4bool isInitialised;

  G4VEnergySpectrum* energySpectrum;
 

  G4PenelopeBremsstrahlungAngular* GetAngularDataForZ(G4int iZ);
  // Map to the objects containing tha angular data
  std::map<G4int,G4PenelopeBremsstrahlungAngular*> *angularData;

  G4PenelopeBremsstrahlungContinuous* GetStoppingPowerData(G4int iZ,G4double energyCut,
							   const G4ParticleDefinition*);
  std::map<std::pair<G4int,G4double>,G4PenelopeBremsstrahlungContinuous*> *stoppingPowerData;


  G4int SampleRandomAtom(const G4MaterialCutsCouple*,G4double energy) const;
  G4VCrossSectionHandler* crossSectionHandler;

 
};

#endif

