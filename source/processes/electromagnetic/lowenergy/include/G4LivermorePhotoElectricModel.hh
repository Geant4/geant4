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
// $Id: G4LivermorePhotoElectricModel.hh,v 1.4 2010-03-15 09:02:29 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Sebastien Incerti
//         30 October 2008
//         on base of G4LowEnergyPhotoElectric developed by A.Forti and M.G.Pia
//
// 15 Mar 2010   L. Pandola, removed methods to set explicitely fluorescence cuts.
//               Main cuts from G4ProductionCutsTable are always used
// 30 May 2011   A Mantero & V Ivanchenko Migration to model design for deexcitation
//


#ifndef G4LivermorePhotoElectricModel_h
#define G4LivermorePhotoElectricModel_h 1

#include "G4VEmModel.hh"

class G4ParticleChangeForGamma;
class G4VCrossSectionHandler;
class G4VAtomDeexcitation;
class G4VPhotoElectricAngularDistribution;
class G4AtomicTransitionManager;

class G4LivermorePhotoElectricModel : public G4VEmModel
{

public:

  G4LivermorePhotoElectricModel(const G4ParticleDefinition* p = 0, 
				const G4String& nam = "LivermorePhElectric");

  virtual ~G4LivermorePhotoElectricModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
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

  //  void ActivateAuger(G4bool);

  void SetAngularGenerator(G4VPhotoElectricAngularDistribution* distribution);

  void SetAngularGenerator(const G4String& name);

protected:

  G4ParticleChangeForGamma* fParticleChange;

private:

  G4ParticleDefinition*     theGamma;
  G4ParticleDefinition*     theElectron;

  G4double lowEnergyLimit;  
  G4double highEnergyLimit; 

  G4int verboseLevel;
  G4bool isInitialised;

  G4VCrossSectionHandler* crossSectionHandler;
  G4VCrossSectionHandler* shellCrossSectionHandler;

  G4VAtomDeexcitation*             fAtomDeexcitation;
  const G4AtomicTransitionManager* fTransitionManager;

  G4VPhotoElectricAngularDistribution* fElectronAngularGenerator;
  G4String generatorName;

  G4LivermorePhotoElectricModel & operator=(const  G4LivermorePhotoElectricModel &right);
  G4LivermorePhotoElectricModel(const  G4LivermorePhotoElectricModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
