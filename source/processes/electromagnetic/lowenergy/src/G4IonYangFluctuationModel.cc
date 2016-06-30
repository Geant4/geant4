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
// 10/05/2001  V.Ivanchenko Clean up againist Linux compilation with -Wall
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

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4Exp.hh"

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

G4double G4IonYangFluctuationModel::TheValue(const G4ParticleDefinition* aParticle,
					     const G4Material* material,
					     G4double kineticEnergy)
{
  G4double mass = aParticle->GetPDGMass() ;
  G4double charge = (aParticle->GetPDGCharge())/eplus ;

  G4double q = YangFluctuationModel(material,kineticEnergy,mass,charge);

  return q ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonYangFluctuationModel::HighEnergyLimit(const G4ParticleDefinition*,
						    const G4Material*) const
{
  return 1.0*TeV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonYangFluctuationModel::LowEnergyLimit(const G4ParticleDefinition* ,
						   const G4Material* ) const
{
  return 0.0 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonYangFluctuationModel::HighEnergyLimit(const G4ParticleDefinition* ) const
{
  return 1.0*TeV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonYangFluctuationModel::LowEnergyLimit(
                              const G4ParticleDefinition* ) const
{
  return 0.0 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4IonYangFluctuationModel::IsInCharge(const G4DynamicParticle*,
					     const G4Material* ) const
{
  return true ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4IonYangFluctuationModel::IsInCharge(const G4ParticleDefinition* ,
					     const G4Material* ) const
{
  return true ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonYangFluctuationModel::YangFluctuationModel(const G4Material* material,
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
  static const G4double b[5][4] = {
  {0.1014,  0.3700,  0.9642,  3.987},
  {0.1955,  0.6941,  2.522,   1.040},
  {0.05058, 0.08975, 0.1419, 10.80},
  {0.05009, 0.08660, 0.2751,  3.787},
  {0.01273, 0.03458, 0.3951,  3.812}
  } ;

  // protons (hadrons)
  if(1.5 > charge) {
    if( kStateGas != material->GetState() ) i = 1 ;

  // ions
  } else {
    G4double zeff = (material->GetElectronDensity())/
                    (material->GetTotNbOfAtomsPerVolume()) ;
    factor = charge * std::pow(charge/zeff, 0.3333) ;

    if( kStateGas == material->GetState() ) {
      energy /= (charge * std::sqrt(charge)) ;

      if(1 == (material->GetNumberOfElements())) {
        i = 2 ;
      } else {
        i = 3 ;
      }

    } else {
      energy /= (charge * std::sqrt(charge*zeff)) ;
      i = 4 ;
    }
  }

  G4double x = b[i][2] * (1.0 - G4Exp( - energy * b[i][3] )) ;

  G4double q = factor * x * b[i][0] /
             ((energy - b[i][1])*(energy - b[i][1]) + x*x) ;

  return q ;
}
