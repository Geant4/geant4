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
// $Id: G4LivermorePolarizedPhotoElectricModel.hh 66241 2012-12-13 18:34:42Z gunter $
//
// Authors: G.Depaola & F.Longo
//

#ifndef G4LivermorePolarizedPhotoElectricModel_h
#define G4LivermorePolarizedPhotoElectricModel_h 1

#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4CrossSectionHandler.hh"

class G4VAtomDeexcitation;

class G4LivermorePolarizedPhotoElectricModel : public G4VEmModel
{

public:

  G4LivermorePolarizedPhotoElectricModel( 
	     const G4String& nam = "LivermorePolarizedPhotoElectric");

  virtual ~G4LivermorePolarizedPhotoElectricModel();

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

protected:

  G4ParticleChangeForGamma* fParticleChange;

private:
  
  G4LivermorePolarizedPhotoElectricModel & operator=(const  G4LivermorePolarizedPhotoElectricModel &right);
  G4LivermorePolarizedPhotoElectricModel(const  G4LivermorePolarizedPhotoElectricModel&);

  G4ParticleDefinition*   theGamma;
  G4ParticleDefinition*   theElectron;

  G4double lowEnergyLimit;  
  G4double highEnergyLimit; 
  G4int    verboseLevel;
  G4bool   fDeexcitationActive;

  G4VCrossSectionHandler* crossSectionHandler;
  G4VCrossSectionHandler* shellCrossSectionHandler;

  G4VAtomDeexcitation*    fAtomDeexcitation;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
