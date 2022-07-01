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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eBremParametrizedModel
//                extention of standard G4eBremsstrahlungModel
//
// Author:        Andreas Schaelicke 
//
// Creation date: 28.03.2008
//
// Modifications:
//
//
// Class Description:
//
// Implementation of energy loss for gamma emission by electrons and
// positrons including an improved version of the LPM effect

// -------------------------------------------------------------------
//

#ifndef G4eBremParametrizedModel_h
#define G4eBremParametrizedModel_h 1

#include "G4VEmModel.hh"
#include "G4NistManager.hh"

class G4ParticleChangeForLoss;
class G4PhysicsVector;

class G4eBremParametrizedModel : public G4VEmModel
{

public:

  explicit G4eBremParametrizedModel(const G4ParticleDefinition* p = nullptr, 
				    const G4String& nam = "eBremParam");

  ~G4eBremParametrizedModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  void InitialiseLocal(const G4ParticleDefinition*, 
		       G4VEmModel* masterModel) override;

  G4double MinEnergyCut(const G4ParticleDefinition*, 
			const G4MaterialCutsCouple*) override;

  G4double ComputeDEDXPerVolume(const G4Material*,
				const G4ParticleDefinition*,
				G4double kineticEnergy,
				G4double cutEnergy) override;
					
  G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
				      G4double tkin, 
				      G4double Z, G4double,
				      G4double cutEnergy,
				      G4double maxEnergy = DBL_MAX) override;
  
  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double cutEnergy,
			 G4double maxEnergy) override;

  void SetupForMaterial(const G4ParticleDefinition*,
                        const G4Material*,G4double) override;

  // hide assignment operator
  G4eBremParametrizedModel & operator=(const  G4eBremParametrizedModel &right) = delete;
  G4eBremParametrizedModel(const  G4eBremParametrizedModel&) = delete;

private:

  void InitialiseConstants();

  G4double ComputeBremLoss(G4double cutEnergy);

  G4double ComputeXSectionPerAtom(G4double cutEnergy);

  G4double ComputeDXSectionPerAtom(G4double gammaEnergy);

  G4double ComputeParametrizedDXSectionPerAtom(G4double kineticEnergy, 
					       G4double gammaEnergy, 
					       G4double Z);

  void SetParticle(const G4ParticleDefinition* p);

  G4double ScreenFunction1(G4double ScreenVariable);

  G4double ScreenFunction2(G4double ScreenVariable);

  inline void SetCurrentElement(const G4double);

protected:

  G4NistManager*              nist;
  const G4ParticleDefinition* particle;
  G4ParticleDefinition*       theGamma;
  G4ParticleChangeForLoss*    fParticleChange;

  static const G4double xgi[8], wgi[8];

  G4double minThreshold;

  // cash
  G4double particleMass;
  G4double kinEnergy;
  G4double totalEnergy;
  G4double currentZ;
  G4double z13, z23, lnZ;
  G4double densityFactor;
  G4double densityCorr;
  G4double Fel, Finel;
  G4double facFel, facFinel;
  G4double fMax,fCoulomb;

private:

  G4double lowKinEnergy;
  G4double fMigdalConstant;
  G4double bremFactor;

  G4bool isInitialised;

protected:

  G4bool isElectron;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4eBremParametrizedModel::SetCurrentElement(const G4double Z)
{
  if(Z != currentZ) {
    currentZ = Z;

    G4int iz = G4lrint(Z);
    z13 = nist->GetZ13(iz);
    z23 = z13*z13;
    lnZ = nist->GetLOGZ(iz);

    Fel = facFel - lnZ/3. ;
    Finel = facFinel - 2.*lnZ/3. ;

    fCoulomb = GetCurrentElement()->GetfCoulomb();
    fMax = Fel-fCoulomb + Finel/currentZ  +  (1.+1./currentZ)/12.;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#endif
