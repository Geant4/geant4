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
// $Id: G4hLDMBremModel.hh 97391 2016-06-02 10:08:45Z gcosmo $
//
// -------------------------------------------------------------------
//
// 21.03.17 V. Grichine based on G4hBremsstrahlungModel
//
// Class Description:
//
// Implementation of energy loss for LDMPhoton emission by hadrons
//
//
//

#ifndef G4LDMBremModel_h
#define G4LDMBremModel_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MuBremsstrahlungModel.hh"

class G4ParticleDefinition;

class G4LDMBremModel : public G4MuBremsstrahlungModel
{

public:

  explicit G4LDMBremModel(const G4ParticleDefinition* p = nullptr,
                         const G4String& nam = "ldmBrem");

  virtual ~G4LDMBremModel();

  inline void SetEpsilon(G4double e){fEpsilon=e;};
  inline G4double GetEpsilon(){return fEpsilon;};

protected:

  virtual G4double ComputeDEDXPerVolume(const G4Material*,
                                const G4ParticleDefinition*,
                                G4double kineticEnergy,
                                G4double cutEnergy) override;

  virtual G4double ComputeDMicroscopicCrossSection( 
                                 G4double tkin,
                                 G4double Z,
                                 G4double gammaEnergy ) override;


  virtual G4double ComputeCrossSectionPerAtom(
                                 const G4ParticleDefinition*,
                                 G4double kineticEnergy,
                                 G4double Z, G4double A,
                                 G4double cutEnergy,
                                 G4double maxEnergy) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin,
                                 G4double maxEnergy) override;

private:

  // hide assignment operator
  G4LDMBremModel & 
    operator=(const  G4LDMBremModel &right) = delete;
  G4LDMBremModel(const  G4LDMBremModel&) = delete;

  const G4ParticleDefinition* theLDMPhoton;

  G4double fEpsilon; // correction for fine_structure_const
  G4double fLDMPhotonMass; 

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
