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
// File name:     G4hNuclearStoppingModel
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 20 July 2000
//
// Modifications: 
// 20/07/2000  V.Ivanchenko First implementation
//
// Class Description: 
//
// Low energy protons/ions ionisation parametrisation
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hNuclearStoppingModel.hh" 
#include "G4UnitsTable.hh"
#include "globals.hh"
#include "G4hZiegler1977Nuclear.hh"
#include "G4hZiegler1985Nuclear.hh"
#include "G4hMollereNuclear.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ElementVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hNuclearStoppingModel::G4hNuclearStoppingModel(const G4String& name):
  G4VLowEnergyModel(name),
  highEnergyLimit(100.*MeV),
  lowEnergyLimit(1.*eV),
  factorPDG2AMU(1.007276/proton_mass_c2),
  theZieglerFactor(eV*cm2*1.0e-15) 
{
  // Registration of parametrisation models of nuclear energy losses
  if("Ziegler1977" == name) { 
      nStopingPowerTable = new G4hZiegler1977Nuclear();
    
  } else if("Ziegler1985" == name) {
      nStopingPowerTable = new G4hZiegler1985Nuclear();

  } else if("ICRU_R49Mollere" == name || " " == name) {
      nStopingPowerTable = new G4hMollereNuclear();
        
  } else {
    G4cout << 
    "G4hLowEnergyIonisation warning: There is no table with the name <" 
 << name << ">" 
 << " for nuclear stopping, <ICRU_R49Mollere> is applied " 
 << G4endl; 
    nStopingPowerTable = new G4hMollereNuclear();
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hNuclearStoppingModel::~G4hNuclearStoppingModel() 
{
  delete nStopingPowerTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hNuclearStoppingModel::TheValue(
                          const G4DynamicParticle* particle,
	       	          const G4Material* material) 
{
  // Projectile nucleus
  G4double energy = particle->GetKineticEnergy() ;
  G4double z1 = abs((particle->GetCharge())/eplus) ;
  G4double m1 = (particle->GetMass())*factorPDG2AMU ;

  G4double nloss = StoppingPower(material, energy, z1, m1) * theZieglerFactor; 

  return nloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hNuclearStoppingModel::TheValue(
                          const G4ParticleDefinition* aParticle,
		          const G4Material* material,
                                G4double kineticEnergy)  
{
  // Projectile nucleus
  G4double z1 = abs((aParticle->GetPDGCharge())/eplus) ;
  G4double m1 = (aParticle->GetPDGMass())*factorPDG2AMU ;

  G4double nloss = StoppingPower(material, kineticEnergy, z1, m1)
                 * theZieglerFactor; 

  return nloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hNuclearStoppingModel::StoppingPower(
                          const G4Material* material,
                                G4double kineticEnergy,
                                G4double z1, G4double m1) const
{
  // Target nucleus
  G4int NumberOfElements = material->GetNumberOfElements() ;
  if(0 == NumberOfElements) return 0.0 ;

  const G4ElementVector* theElementVector = 
                                 material->GetElementVector() ;
  const G4double* theAtomicNumDensityVector = 
                                 material->GetAtomicNumDensityVector() ;

  //  loop for the elements in the material

  G4double nloss = 0.0;

  for (G4int iel=0; iel<NumberOfElements; iel++) {
    const G4Element* element = (*theElementVector)(iel) ;
    G4double z2 = element->GetZ() ;
    G4double m2 = element->GetA()*mole/g ;
    nloss += (nStopingPowerTable->
              NuclearStoppingPower(z1, z2, m1, m2, kineticEnergy))
           * theAtomicNumDensityVector[iel] ;    
  }

  return nloss;
}









