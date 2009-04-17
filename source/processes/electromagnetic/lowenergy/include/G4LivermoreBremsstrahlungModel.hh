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
// $Id: G4LivermoreBremsstrahlungModel.hh,v 1.2 2009-04-17 10:29:20 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 03 Mar 2009   L. Pandola   1st implementation. Migration from EM process 
//                            to EM model. Physics is unchanged.
//
// -------------------------------------------------------------------
//
// Class description:
// Low Energy Electromagnetic Physics, e- bremsstrahlung
// with Livermore Model
// -------------------------------------------------------------------

#ifndef G4LIVERMOREBREMSSTRAHLUNGMODEL_HH
#define G4LIVERMOREBREMSSTRAHLUNGMODEL_HH 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4ParticleChangeForLoss.hh"

class G4ParticleDefinition;
class G4MaterialCutsCouple;
class G4Material;
class G4VBremAngularDistribution;
class G4BremsstrahlungCrossSectionHandler;
class G4VEnergySpectrum;

class G4LivermoreBremsstrahlungModel : public G4VEmModel 
{

public:
  
  G4LivermoreBremsstrahlungModel(const G4ParticleDefinition* p=0,
			 const G4String& processName = "LowEnBrem");
  
  virtual ~G4LivermoreBremsstrahlungModel();

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

  void SetVerboseLevel(G4int vl) {verboseLevel = vl;};

  void SetAngularGenerator(G4VBremAngularDistribution* distribution);
  void SetAngularGenerator(const G4String& name);

protected:
  G4ParticleChangeForLoss* fParticleChange;

private:
 
  G4LivermoreBremsstrahlungModel & operator=(const G4LivermoreBremsstrahlungModel &right);
  G4LivermoreBremsstrahlungModel(const G4LivermoreBremsstrahlungModel&);


  //Intrinsic energy limits of the model: cannot be extended by the parent process
  G4double fIntrinsicLowEnergyLimit;
  G4double fIntrinsicHighEnergyLimit;
  G4int fNBinEnergyLoss;

  G4bool isInitialised;
 
  G4int verboseLevel;

  G4BremsstrahlungCrossSectionHandler* crossSectionHandler;
  G4VEnergySpectrum* energySpectrum;
  G4DataVector energyBins;

  G4VBremAngularDistribution* angularDistribution;
  G4VBremAngularDistribution* TsaiAngularDistribution;
  G4String generatorName;
};

#endif

