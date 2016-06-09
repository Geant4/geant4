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
// $Id: G4PEEffectModel.hh,v 1.6 2007/05/22 17:34:36 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4PEEffectModel
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 21.04.2005
//
// Modifications:
//
// 06.02.2006 : added ComputeMeanFreePath()  (mma)
//
// Class Description:
//
// Implementation of the photo-electric effect
//

// -------------------------------------------------------------------
//

#ifndef G4PEEffectModel_h
#define G4PEEffectModel_h 1

#include "G4VEmModel.hh"
#include "G4PhysicsTable.hh"

class G4ParticleChangeForGamma;

class G4PEEffectModel : public G4VEmModel
{

public:

  G4PEEffectModel(const G4ParticleDefinition* p = 0,
		  const G4String& nam = "PhotoElectric");

  virtual ~G4PEEffectModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                      G4double kinEnergy,
                                      G4double Z,
                                      G4double A,
                                      G4double, G4double);
				      
  G4double ComputeMeanFreePath( const G4ParticleDefinition*,
                                      G4double kinEnergy,
				const G4Material* material,      
                                      G4double, G4double);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);
protected:

  virtual G4double ElecCosThetaDistribution(G4double ElecKineEnergy);

private:

  G4ParticleDefinition*     theGamma;
  G4ParticleDefinition*     theElectron;
  G4ParticleChangeForGamma* fParticleChange;

  G4double                  fminimalEnergy;
  G4bool                    isInitialized;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

inline G4double G4PEEffectModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double energy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
 G4double* SandiaCof = G4SandiaTable::GetSandiaCofPerAtom((G4int)Z, energy);

 G4double energy2 = energy*energy, energy3 = energy*energy2,
          energy4 = energy2*energy2;

 return SandiaCof[0]/energy  + SandiaCof[1]/energy2 +
        SandiaCof[2]/energy3 + SandiaCof[3]/energy4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4PEEffectModel::ComputeMeanFreePath(
                                       const G4ParticleDefinition*,
                                             G4double energy,
				       const G4Material* material,
                                             G4double, G4double)
{
 G4double* SandiaCof = material->GetSandiaTable()
                                ->GetSandiaCofForMaterial(energy);
				
 G4double energy2 = energy*energy, energy3 = energy*energy2,
          energy4 = energy2*energy2;
	  
 G4double cross = SandiaCof[0]/energy  + SandiaCof[1]/energy2 +
                  SandiaCof[2]/energy3 + SandiaCof[3]/energy4; 
 
 G4double mfp = DBL_MAX;
 if (cross > 0.) mfp = 1./cross;
 return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
