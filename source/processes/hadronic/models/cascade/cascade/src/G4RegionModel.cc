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
//
// 20100319  M. Kelsey -- Eliminate unnecessary use of std::pow()
// 20101019  M. Kelsey -- CoVerity report: unitialized constructor

#include "G4RegionModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4HadronicException.hh"
#include "G4InuclSpecialFunctions.hh"

using namespace G4InuclSpecialFunctions;

const G4double G4RegionModel::radius0 = 1.0E-15; 
const G4double G4RegionModel::BE = 7;

G4RegionModel::G4RegionModel(const G4int numberOfLayers,
			     const G4int A, const G4int Z)
  : massNumber(A), protonNumber(Z)
{
  //count the radiuses, densities and fermi momenta with A and Z
  G4double r = radius0*G4cbrt(A);

  if(numberOfLayers==1){ 
    radius.push_back(r);

    G4double vol = 4.0/3.0 * pi * r*r*r;
    G4double rho = G4double(A) / vol;
    density.push_back(rho);

    G4double protonMass = G4Proton::Proton()->GetPDGMass();
    G4double neutronMass = G4Neutron::Neutron()->GetPDGMass();
    G4double protonDensity = G4double(Z) / vol;
    G4double neutronDensity = G4double(A-Z) / vol;

    protonFermiEnergy.push_back(GetFermiEnergy(protonDensity, protonMass));
    neutronFermiEnergy.push_back(GetFermiEnergy(neutronDensity, neutronMass));
    
    protonFermiMomentum.push_back(GetFermiMomentum(protonDensity, protonMass));
    neutronFermiMomentum.push_back(GetFermiMomentum(neutronDensity, neutronMass));

    G4double fermiEP = *protonFermiEnergy.begin();
    G4double fermiEN = *neutronFermiEnergy.begin();
    protonPotentialEnergy.push_back(-(fermiEP + BE));
    neutronPotentialEnergy.push_back(-(fermiEN + BE));
  }
  else{
  if(numberOfLayers==3){
    radius.push_back(0.1*r);
    radius.push_back(0.2*r);
    radius.push_back(0.9*r);
    
  }
  }
}

G4RegionModel::~G4RegionModel(){}

G4double G4RegionModel::GetDensity(G4double r){
  my_iterator j=density.begin();
     for(my_iterator i=radius.begin(); i<radius.end(); i++){
     if(r <= *i) return *j;
     j++;
   }
   return 0;
}

G4double G4RegionModel::GetPotentialEnergy(G4double r, G4int particle){
  if(particle == 0){ //proton
    my_iterator j=protonPotentialEnergy.begin();
     for(my_iterator i=radius.begin(); i<radius.end(); i++){
     if(r <= *i) return *j;
     j++;
   }
     return 0;
  }
  
  if(particle == 1){ //neutron
    my_iterator j=neutronPotentialEnergy.begin();
     for(my_iterator i=radius.begin(); i<radius.end(); i++){
     if(r <= *i) return *j;
     j++;
   }
     return 0;
  }
  return 0;
}

G4double G4RegionModel::GetMaximumNucleonMomentum(G4double r,
						  G4int nucleon){ 
  if(nucleon == 0){
     my_iterator j=protonFermiMomentum.begin();
     for(my_iterator i=radius.begin(); i<radius.end(); i++){
     if(r <= *i)  return *j;
     j++;
     }
  }
  if(nucleon==1){
     my_iterator j=neutronFermiMomentum.begin();
     for(my_iterator i=radius.begin(); i<radius.end(); i++){
     if(r <= *i)  return *j;
     j++;
     }
  }
  throw G4HadronicException(__FILE__, __LINE__, "G4RegionModel::GetMaximumNucleonMomentum - return value undefined");
  return 0;

}

G4double G4RegionModel::GetFermiMomentum(G4double aDensity,
					 G4double aMass){
  return std::sqrt(2*aMass*GetFermiEnergy(aDensity, aMass));
}

G4double G4RegionModel::GetFermiEnergy(G4double aDensity,
					 G4double aMass){
  G4double densFactor = G4cbrt(3.0*pi2*aDensity);		// 2/3 power
  densFactor *= densFactor;

  return hbar_Planck*hbar_Planck/(2.0*aMass) * densFactor;
}














