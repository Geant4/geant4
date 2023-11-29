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
// GEANT4 Class file
//
//
// File name:     G4hNuclearStoppingModel
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 20 July 2000
//
// Modifications:
// 20/07/2000  V.Ivanchenko First implementation
// 22/08/2000  V.Ivanchenko Bug fixed in call of a model
// 03/10/2000  V.Ivanchenko CodeWizard clean up
//
// Class Description: 
//
// Low energy protons/ions nuclear stopping parametrisation
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hNuclearStoppingModel.hh" 

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4hZiegler1985Nuclear.hh"
#include "G4hICRU49Nuclear.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ElementVector.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hNuclearStoppingModel::G4hNuclearStoppingModel(const G4String& name)
  :G4VLowEnergyModel(name), modelName(name)
{
  InitializeMe() ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hNuclearStoppingModel::InitializeMe()
{
  // Constants
  highEnergyLimit = 100.*MeV ;
  lowEnergyLimit  = 1.*eV ;
  factorPDG2AMU   = 1.007276/proton_mass_c2 ;
  theZieglerFactor= eV*cm2*1.0e-15 ; 

  // Registration of parametrisation models of nuclear energy losses
  G4String blank  = G4String(" ") ;
  G4String ir49   = G4String("ICRU_R49") ;
  G4String zi85   = G4String("Ziegler1985") ;
  if(ir49 == modelName || blank == modelName) {
      nStopingPowerTable = new G4hICRU49Nuclear();

  } else if(zi85 == modelName) {
      nStopingPowerTable = new G4hZiegler1985Nuclear();
        
  } else {
    G4cout << 
    "G4hLowEnergyIonisation warning: There is no table with the modelName <" 
	   << modelName << ">" 
	   << " for nuclear stopping, <ICRU_R49> is applied " 
	   << G4endl; 
    nStopingPowerTable = new G4hICRU49Nuclear();
  }

  // Default is nuclear stopping fluctuations On
  //  nStopingPowerTable->SetNuclearStoppingFluctuationsOn();  
  nStopingPowerTable->SetNuclearStoppingFluctuationsOff();  
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
  G4double z1 = std::abs((particle->GetCharge())/eplus) ;
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
  G4double z1 = std::abs((aParticle->GetPDGCharge())/eplus) ;
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
  std::size_t NumberOfElements = material->GetNumberOfElements() ;
  if(0 == NumberOfElements) return 0.0 ;

  const G4ElementVector* theElementVector = 
                                 material->GetElementVector() ;
  const G4double* theAtomicNumDensityVector = 
                                 material->GetAtomicNumDensityVector() ;

  //  loop for the elements in the material

  G4double nloss = 0.0;

  for (std::size_t iel=0; iel<NumberOfElements; ++iel) {
    const G4Element* element = (*theElementVector)[iel] ;
    G4double z2 = element->GetZ();
    G4double m2Local = element->GetA()*mole/g ;
    nloss += (nStopingPowerTable->
              NuclearStoppingPower(kineticEnergy, z1, z2, m1, m2Local))
           * theAtomicNumDensityVector[iel] ;    
  }

  return nloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
