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
// $Id$
// GEANT4 tag $Name: geant4-09-04 $
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

  G4eBremParametrizedModel(const G4ParticleDefinition* p = 0, 
			   const G4String& nam = "eBremParam");

  virtual ~G4eBremParametrizedModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double MinEnergyCut(const G4ParticleDefinition*, 
				const G4MaterialCutsCouple*);

  virtual G4double ComputeDEDXPerVolume(const G4Material*,
					const G4ParticleDefinition*,
					G4double kineticEnergy,
					G4double cutEnergy);
					
  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
					      G4double tkin, 
					      G4double Z,   G4double,
					      G4double cutEnergy,
					      G4double maxEnergy = DBL_MAX);
  
  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double cutEnergy,
				 G4double maxEnergy);

  virtual void SetupForMaterial(const G4ParticleDefinition*,
                                const G4Material*,G4double);


private:

  void InitialiseConstants();

  G4double ComputeBremLoss(G4double cutEnergy);

  G4double ComputeXSectionPerAtom(G4double cutEnergy);

  G4double ComputeDXSectionPerAtom(G4double gammaEnergy);

  void SetParticle(const G4ParticleDefinition* p);

  // * fast inline functions *
  inline void SetCurrentElement(const G4double);

  // hide assignment operator
  G4eBremParametrizedModel & operator=(const  G4eBremParametrizedModel &right);
  G4eBremParametrizedModel(const  G4eBremParametrizedModel&);

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

  G4bool   isElectron;

private:

  // consts
  G4double lowKinEnergy;
  G4double fMigdalConstant;
  G4double bremFactor;

  G4bool   isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4eBremParametrizedModel::SetCurrentElement(const G4double Z)
{
  std::cout<<"SetCurrentElement Z="<<Z<<std::endl;
  if(Z != currentZ) {
    currentZ = Z;

    G4int iz = G4int(Z);
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
