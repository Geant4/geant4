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
// $Id: G4teoCrossSection.cc,v 1.9 2011-01-03 19:35:11 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//         
//
// History:
// -----------
//  21 Apr 2009   ALF   1st implementation
//  29 Apr 2009   ALF Updated Desing for Integration
//


#include "globals.hh"
#include "G4teoCrossSection.hh"
#include "G4Proton.hh"
//#include "G4Alpha.hh"

G4teoCrossSection::G4teoCrossSection(const G4String& nam)
  :G4VhShellCrossSection(nam),totalCS(0.0)
{ 
  ecpssrShellK  = new G4AnalyticalEcpssrKCrossSection();  
  ecpssrShellLi = new G4AnalyticalEcpssrLiCrossSection();
}

G4teoCrossSection::~G4teoCrossSection()
{ 
  delete ecpssrShellK;
  delete ecpssrShellLi;
}

std::vector<G4double> G4teoCrossSection::GetCrossSection(G4int Z,
							 G4double incidentEnergy,
							 G4double mass,
							 G4double,
							 G4bool) const
{
  std::vector<G4double> crossSections;

  crossSections.push_back( ecpssrShellK->CalculateCrossSection(Z, mass, incidentEnergy) );
  
  crossSections.push_back( ecpssrShellLi->CalculateL1CrossSection(Z, mass, incidentEnergy) );
  crossSections.push_back( ecpssrShellLi->CalculateL2CrossSection(Z, mass, incidentEnergy) );
  crossSections.push_back( ecpssrShellLi->CalculateL3CrossSection(Z, mass, incidentEnergy) );

  return crossSections;
}

G4double G4teoCrossSection::CrossSection(G4int Z, G4int shell,
					 G4double incidentEnergy,
					 G4double mass) const
{
  G4double res = 0.0;
  if(3 < shell) {
    return res; 
  } else if(0 == shell) { 
    res = ecpssrShellK->CalculateCrossSection(Z, mass, incidentEnergy);
  } else if(1 == shell) { 
    res = ecpssrShellLi->CalculateL1CrossSection(Z, mass, incidentEnergy);
  } else if(2 == shell) { 
    res = ecpssrShellLi->CalculateL2CrossSection(Z, mass, incidentEnergy);
  } else if(3 == shell) { 
    res = ecpssrShellLi->CalculateL3CrossSection(Z, mass, incidentEnergy);
  }
  return res;
}

std::vector<G4double> G4teoCrossSection::Probabilities(G4int Z,
						       G4double incidentEnergy,
						       G4double mass,
						       G4double deltaEnergy) const
{
  std::vector<G4double> crossSections = 
    GetCrossSection(Z, incidentEnergy, mass, deltaEnergy);

  for (size_t i=0; i<crossSections.size(); i++ ) {
    
    if (totalCS) {
      crossSections[i] = crossSections[i]/totalCS;
    }
    
  }
  return crossSections;
}

void G4teoCrossSection::SetTotalCS(G4double val){

  totalCS = val;
  //  G4cout << "totalXS set to: " << val / barn << " barns" << G4endl;
}



