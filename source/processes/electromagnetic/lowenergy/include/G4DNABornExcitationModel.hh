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
// $Id: G4DNABornExcitationModel.hh,v 1.1.4.1 2010/04/01 09:07:24 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-03-patch-02 $
//

#ifndef G4DNABornExcitationModel_h
#define G4DNABornExcitationModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"

#include "G4DNACrossSectionDataSet.hh"
#include "G4LogLogInterpolation.hh"
#include "G4Proton.hh"
#include "G4WaterExcitationStructure.hh"

class G4DNABornExcitationModel : public G4VEmModel
{

public:

  G4DNABornExcitationModel(const G4ParticleDefinition* p = 0, 
		          const G4String& nam = "DNABornExcitationModel");

  virtual ~G4DNABornExcitationModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double CrossSectionPerVolume(const G4Material* material,
					   const G4ParticleDefinition* p,
					   G4double ekin,
					   G4double emin,
					   G4double emax);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:

  G4double lowEnergyLimit;
  G4double highEnergyLimit;
  G4bool isInitialised;
  G4int verboseLevel;
  
  // Cross section

  G4DNACrossSectionDataSet* table;

  // Partial cross section
  
  G4int RandomSelect(G4double energy);
  
  // Final state

  G4WaterExcitationStructure waterStructure;
   
  //
  
  G4DNABornExcitationModel & operator=(const  G4DNABornExcitationModel &right);
  G4DNABornExcitationModel(const  G4DNABornExcitationModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
