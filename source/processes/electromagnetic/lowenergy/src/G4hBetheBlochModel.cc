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
// File name:     G4hBetheBlochModel
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 20 July 2000
//
// Modifications: 
// 20/07/2000  V.Ivanchenko First implementation
// 03/10/2000  V.Ivanchenko clean up accoding to CodeWizard
//
// Class Description: 
//
// Bethe-Bloch ionisation model
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hBetheBlochModel.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hBetheBlochModel::G4hBetheBlochModel(const G4String& name)
  : G4VLowEnergyModel(name), 
    lowEnergyLimit(1.*MeV),
    highEnergyLimit(100.*GeV),
    twoln10(2.*std::log(10.)),
    bg2lim(0.0169), 
    taulim(8.4146e-3)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hBetheBlochModel::~G4hBetheBlochModel() 
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hBetheBlochModel::TheValue(const G4DynamicParticle* particle,
	       	                      const G4Material* material) 
{
  G4double energy = particle->GetKineticEnergy() ;
  G4double particleMass = particle->GetMass() ;

  G4double eloss  = BetheBlochFormula(material,energy,particleMass) ;

  return eloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hBetheBlochModel::TheValue(const G4ParticleDefinition* aParticle,
       		                      const G4Material* material,
				      G4double kineticEnergy) 
{
  G4double particleMass = aParticle->GetPDGMass() ;
  G4double eloss  = BetheBlochFormula(material,kineticEnergy,particleMass) ;

  return eloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hBetheBlochModel::HighEnergyLimit(
                             const G4ParticleDefinition* ,
                             const G4Material* ) const
{
  return highEnergyLimit ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hBetheBlochModel::LowEnergyLimit(
                             const G4ParticleDefinition* aParticle,
                             const G4Material* material) const
{
  G4double taul = (material->GetIonisation()->GetTaul())*
                  (aParticle->GetPDGMass()) ;
  return taul ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hBetheBlochModel::HighEnergyLimit(
                             const G4ParticleDefinition* ) const
{
  return highEnergyLimit ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hBetheBlochModel::LowEnergyLimit(
                             const G4ParticleDefinition* ) const
{
  return lowEnergyLimit ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4bool G4hBetheBlochModel::IsInCharge(const G4DynamicParticle* ,
		                      const G4Material* ) const
{
  return true ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4bool G4hBetheBlochModel::IsInCharge(const G4ParticleDefinition* ,
		                      const G4Material* ) const
{
  return true ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hBetheBlochModel::BetheBlochFormula(
                             const G4Material* material,
                                   G4double kineticEnergy,
                                   G4double particleMass) const
{
  // This member function is applied normally to proton/antiproton
  G4double ionloss ;

  G4double rateMass = electron_mass_c2/particleMass ;

  G4double taul = material->GetIonisation()->GetTaul() ;
  G4double tau  = kineticEnergy/particleMass ;    // tau is relative energy
  
  // It is not normal case for this function
  // for low energy parametrisation have to be applied
  if ( tau < taul ) tau = taul ; 
    
  // some local variables 
    
  G4double gamma,bg2,beta2,tmax,x,delta,sh ;
  G4double electronDensity = material->GetElectronDensity();
  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double eexc2 = eexc*eexc ;
  G4double cden  = material->GetIonisation()->GetCdensity();
  G4double mden  = material->GetIonisation()->GetMdensity();
  G4double aden  = material->GetIonisation()->GetAdensity();
  G4double x0den = material->GetIonisation()->GetX0density();
  G4double x1den = material->GetIonisation()->GetX1density();
  G4double* shellCorrectionVector =
            material->GetIonisation()->GetShellCorrectionVector();
    
  gamma = tau + 1.0 ;
  bg2 = tau*(tau+2.0) ;
  beta2 = bg2/(gamma*gamma) ;
  tmax = 2.*electron_mass_c2*bg2/(1.+2.*gamma*rateMass+rateMass*rateMass) ;
        
  ionloss = std::log(2.0*electron_mass_c2*bg2*tmax/eexc2)-2.0*beta2 ;
    
  // density correction     
  x = std::log(bg2)/twoln10 ;
  if ( x < x0den ) {
    delta = 0.0 ;

  } else { 
    delta = twoln10*x - cden ;
    if ( x < x1den ) delta += aden*std::pow((x1den-x),mden) ;
  }
    
  // shell correction 
  sh = 0.0 ;      
  x  = 1.0 ;

  if ( bg2 > bg2lim ) {
    for (G4int k=0; k<=2; k++) {
	x *= bg2 ;
	sh += shellCorrectionVector[k]/x;
    }

  } else {
    for (G4int k=0; k<=2; k++) {
	x *= bg2lim ;
	sh += shellCorrectionVector[k]/x;
    }
    sh *= std::log(tau/taul)/std::log(taulim/taul) ;     
  }
    
  // now compute the total ionization loss
    
  ionloss -= delta + sh ;
  ionloss *= twopi_mc2_rcl2*electronDensity/beta2 ;

  if ( ionloss < 0.0) ionloss = 0.0 ;
  
  return ionloss;
}


