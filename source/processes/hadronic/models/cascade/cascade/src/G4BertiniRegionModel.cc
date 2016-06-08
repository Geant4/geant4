//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************

#include "G4BertiniRegionModel.hh"

const G4double G4BertiniRegionModel::radius0 = 1.0E-15; 
const G4double G4BertiniRegionModel::BE = 7;

G4BertiniRegionModel::G4BertiniRegionModel(const G4int numberOfLayers, const G4int A, const G4int Z)
{
  //count the radiuses, densities and fermi momenta with A and Z
  G4double r = radius0*pow(A, 1/3);

  if(numberOfLayers==1){ 
    radius.push_back(r);

    G4double rho = A / (4/3*pi*pow(r,3));
    density.push_back(rho);

    G4double protonMass = G4Proton::Proton()->GetPDGMass();
    G4double neutronMass = G4Neutron::Neutron()->GetPDGMass();
    G4double protonDensity = Z / (4/3*pi*pow(r,3));
    G4double neutronDensity = (A-Z) / (4/3*pi*pow(r,3));

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
G4BertiniRegionModel::~G4BertiniRegionModel(){}

/*
void G4BertiniRegionModel::Init(const G4int numberOfLayers, const G4int A, const G4int Z){
  //count the radiuses, densities and fermi momenta with A and Z
  G4double r = radius0*pow(A, 1/3);

  if(numberOfLayers==1){ 
    radius.push_back(r);

    G4double rho = A / (4/3*pi*pow(r,3)); 
    density.push_back(rho);

    G4double protonMass = G4Proton::Proton()->GetPDGMass()/MeV;
    G4double neutronMass = G4Neutron::Neutron()->GetPDGMass()/MeV ;
    G4double protonDensity = Z / (4/3*pi*pow(r,3));
    G4double neutronDensity = (A-Z) / (4/3*pi*pow(r,3));

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
*/

//--------------------------------------------------------------

G4double G4BertiniRegionModel::GetDensity(G4double r){
  iterator j=density.begin();
     for(iterator i=radius.begin(); i<radius.end(); i++){
     if(r <= *i) return *j;
     j++;
   }
   return 0;
}


G4double G4BertiniRegionModel::GetPotentialEnergy(G4double r, G4int particle){
  if(particle == 0){ //proton
    iterator j=protonPotentialEnergy.begin();
     for(iterator i=radius.begin(); i<radius.end(); i++){
     if(r <= *i) return *j;
     j++;
   }
     return 0;
  }
  
  if(particle == 1){ //neutron
    iterator j=neutronPotentialEnergy.begin();
     for(iterator i=radius.begin(); i<radius.end(); i++){
     if(r <= *i) return *j;
     j++;
   }
     return 0;
  }
  return 0;
}

G4double G4BertiniRegionModel::GetMaximumNucleonMomentum(G4double r,
						  G4int nucleon){
 
  if(nucleon == 0){
     iterator j=protonFermiMomentum.begin();
     for(iterator i=radius.begin(); i<radius.end(); i++){
     if(r <= *i)  return *j;
     j++;
     }
  }
  if(nucleon==1){
     iterator j=neutronFermiMomentum.begin();
     for(iterator i=radius.begin(); i<radius.end(); i++){
     if(r <= *i)  return *j;
     j++;
     }
  }
  G4Exception("G4BertiniRegionModel::GetMaximumNucleonMomentum - return value undefined");
  return 0;

}


G4double G4BertiniRegionModel::GetFermiMomentum(G4double aDensity,
					 G4double aMass){
  
  return sqrt(2*aMass*GetFermiEnergy(aDensity, aMass));
  
 
}

G4double G4BertiniRegionModel::GetFermiEnergy(G4double aDensity,
					 G4double aMass){
  
  //G4double hbar = 1.0E-6;
    return (pow(hbar_Planck,2)/(2*aMass)*pow((3*pi2*aDensity),2/3)); 
 
 
}














