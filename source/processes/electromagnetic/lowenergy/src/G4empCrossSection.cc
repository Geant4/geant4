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
// $Id: G4empCrossSection.cc 92765 2015-09-15 15:20:47Z gcosmo $
//
//         
//
// History:
// -----------
//  29 Apr 2009   ALF   1st implementation
//  15 Mar 2011   ALF introduced the usage of G4AtomicShellEnumerator
//  09 Mar 2012   LP  updated methods
//


#include "globals.hh"
#include "G4empCrossSection.hh"
#include "G4Proton.hh"


G4empCrossSection::G4empCrossSection(const G4String& nam)
  :G4VhShellCrossSection(nam),totalCS(0.0)
{ 
 if (nam == "Empirical") 
  {
    paulShellK = new G4PaulKxsModel();
    orlicShellLi = new G4OrlicLiXsModel();
    flag=0;
  }
  else 
  { 
    G4cout << "G4empCrossSection::G4empCrossSection: " 
	   << "ERROR in G4empCrossSection name; Paul+Orlic is selected." 
	   << G4endl; 
    paulShellK = new G4PaulKxsModel();
    orlicShellLi = new G4OrlicLiXsModel();
    flag=0;
  }

}

G4empCrossSection::~G4empCrossSection()
{ 
  delete paulShellK;
  delete orlicShellLi;
}

std::vector<G4double> G4empCrossSection::GetCrossSection(G4int Z,
							 G4double incidentEnergy,
							 G4double mass,
							 G4double,
							 const G4Material*) 
{
  std::vector<G4double> crossSections;
  G4ParticleDefinition* aProton = G4Proton::Proton();

  crossSections.push_back( paulShellK->CalculateKCrossSection(Z, mass, incidentEnergy) );

  // this check should be done in the Orlic class, that can handle only protons;
  // however this would lead up tp three checks of the mass, while here we have only one
  // moreover, at the present time,this class handles explicitly Paul and Orlic models,
  // so it can hadle the responsibility of this check too

  if (mass == aProton->GetPDGMass()) {

    if (flag==0)
    {
      crossSections.push_back( orlicShellLi->CalculateL1CrossSection(Z, incidentEnergy) );
      crossSections.push_back( orlicShellLi->CalculateL2CrossSection(Z, incidentEnergy) );
      crossSections.push_back( orlicShellLi->CalculateL3CrossSection(Z, incidentEnergy) );
    }

  }

  else {
    crossSections.push_back( 0. );
    crossSections.push_back( 0. );
    crossSections.push_back( 0. );
  }  
  return crossSections;

}

G4double G4empCrossSection::CrossSection(G4int Z, G4AtomicShellEnumerator shell,
					 G4double incidentEnergy,
					 G4double mass,
					 const G4Material*)
{
  G4double res = 0.0;
  G4ParticleDefinition* aProton = G4Proton::Proton();
  
  if(fKShell == shell) { 
    res = paulShellK->CalculateKCrossSection(Z, mass, incidentEnergy);
  } 
  // this check should be done in the Orlic class, that can handle only protons;
  // however this would lead up tp three checks of the mass, while here we have only one
  // moreover, at the present time,this class handles explicitly Paul and Orlic models,
  // so it can hadle the responsibility of this check too

  else if (mass == aProton->GetPDGMass()) {
    
    if(fL1Shell == shell) { 
      if (flag==0) res =   orlicShellLi->CalculateL1CrossSection(Z, incidentEnergy);
    } 
    else if(fL2Shell == shell) { 
      if (flag==0) res =   orlicShellLi->CalculateL2CrossSection(Z, incidentEnergy);
    } 
    else if(fL3Shell == shell) { 
      if (flag==0) res =   orlicShellLi->CalculateL3CrossSection(Z, incidentEnergy);
    } 
  }
  return res;
}

std::vector<G4double> G4empCrossSection::Probabilities(G4int Z,
						       G4double incidentEnergy,
						       G4double mass,
						       G4double deltaEnergy,
						       const G4Material* mat)
{
  
  std::vector<G4double> crossSections = GetCrossSection(Z, incidentEnergy, mass, deltaEnergy,mat);

  for (size_t i=0; i<crossSections.size(); i++ ) {
    
    if (totalCS) {
      crossSections[i] = crossSections[i]/totalCS;
    }
    
  }

  return crossSections;

}


void G4empCrossSection::SetTotalCS(G4double val){

  totalCS = val;

}



