// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// For information related to this code contact:
// Geant4 Collaboration
//
// File name:     G4IonYangFluctuationModel
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 18 August 2000
//
// Modifications: 
// 18/08/2000  V.Ivanchenko First implementation
// 04/09/2000  V.Ivanchenko Rename fluctuations            
// 03/10/2000  V.Ivanchenko CodeWizard clean up
//
// -------------------------------------------------------------------
// Class Description: 
//
// The aproximation of additional ion energy loss fluctuations 
// Q.Yang et al., NIM B61(1991)149-155.
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4IonYangFluctuationModel.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4IonYangFluctuationModel::G4IonYangFluctuationModel(const G4String& name)
  : G4VLowEnergyModel(name)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4IonYangFluctuationModel::~G4IonYangFluctuationModel() 
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonYangFluctuationModel::TheValue(const G4DynamicParticle* particle,
                	                    const G4Material* material) 
{
  G4double energy = particle->GetKineticEnergy() ;
  G4double mass = particle->GetMass() ;
  G4double charge = (particle->GetCharge())/eplus ;

  G4double q = YangFluctuationModel(material,energy,mass,charge) ;

  return q ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonYangFluctuationModel::TheValue(
                                   const G4ParticleDefinition* aParticle,
       		                   const G4Material* material,
                                         G4double kineticEnergy) 
{
  G4double mass = aParticle->GetPDGMass() ;
  G4double charge = (aParticle->GetPDGCharge())/eplus ;

  G4double q = YangFluctuationModel(material,kineticEnergy,mass,charge);

  return q ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonYangFluctuationModel::HighEnergyLimit(
                             const G4ParticleDefinition* aParticle,
                             const G4Material* material) const
{
  return 1.0*TeV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonYangFluctuationModel::LowEnergyLimit(
                              const G4ParticleDefinition* aParticle,
                              const G4Material* material) const
{
  return 0.0 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonYangFluctuationModel::HighEnergyLimit(
                              const G4ParticleDefinition* aParticle) const
{
  return 1.0*TeV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonYangFluctuationModel::LowEnergyLimit(
                              const G4ParticleDefinition* aParticle) const
{
  return 0.0 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4bool G4IonYangFluctuationModel::IsInCharge(
                                  const G4DynamicParticle* particle,
    		                  const G4Material* material) const
{
  return true ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4bool G4IonYangFluctuationModel::IsInCharge(
                                  const G4ParticleDefinition* aParticle,
      		                  const G4Material* material) const
{
  return true ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonYangFluctuationModel::YangFluctuationModel(
                                    const G4Material* material,
                                          G4double kineticEnergy,
                                          G4double particleMass,
                                          G4double charge) const
{
  // The aproximation of energy loss fluctuations 
  // Q.Yang et al., NIM B61(1991)149-155.

  // Reduced energy in MeV/AMU
  G4double energy = kineticEnergy *amu_c2/(particleMass*MeV) ;

  G4int i = 0 ;
  G4double factor = 1.0 ;

  // The index of set of parameters i = 0 for protons(hadrons) in gases
  //                                    1 for protons(hadrons) in solids
  //                                    2 for ions in atomic gases
  //                                    3 for ions in molecular gases
  //                                    4 for ions in solids
  static G4double b[5][4] = {
  0.1014,  0.3700,  0.9642,  3.987,
  0.1955,  0.6941,  2.522,   1.040,
  0.05058, 0.08975, 0.1419, 10.80,
  0.05009, 0.08660, 0.2751,  3.787,
  0.01273, 0.03458, 0.3951,  3.812
  } ;

  // protons (hadrons) 
  if(1.5 > charge) {
    if( kStateGas != material->GetState() ) i = 1 ;

  // ions
  } else {
    G4double zeff = (material->GetElectronDensity())/
                    (material->GetTotNbOfAtomsPerVolume()) ;
    factor = charge * pow(charge/zeff, 0.3333) ;

    if( kStateGas == material->GetState() ) {
      energy /= (charge * sqrt(charge)) ;      

      if(1 == (material->GetNumberOfElements())) {
        i = 2 ;
      } else {
        i = 3 ;
      }

    } else {
      energy /= (charge * sqrt(charge*zeff)) ;      
      i = 4 ;
    }
  }

  G4double x = b[i][2] * (1.0 - exp( - energy * b[i][3] )) ;

  G4double q = factor * x * b[i][0] / 
             ((energy - b[i][1])*(energy - b[i][1]) + x*x) ;

  return q ;    
}






