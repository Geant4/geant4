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
//$Id: G4empCrossSection.cc,v 1.3 2010/11/12 18:09:44 mantero Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//         
//
// History:
// -----------
//  29 Apr 2009   ALF   1st implementation
//
// -------------------------------------------------------------------
// Class description:
// empirical model for K and L Ionization CS for Protons and Alpha
// Further documentation available from http://www.ge.infn.it/geant4/lowE
// -------------------------------------------------------------------


#include "globals.hh"
#include "G4empCrossSection.hh"
#include "G4Proton.hh"
//#include "G4Alpha.hh"
//#include <math.h>

G4empCrossSection::G4empCrossSection()
  :totalCS(0)
{ 

  paulShellK = new G4PaulKCrossSection();
  orlicShellLi = new G4OrlicLiCrossSection();

}

G4empCrossSection::~G4empCrossSection()
{ 

  delete paulShellK;
  delete orlicShellLi;

}

std::vector<G4double> G4empCrossSection::GetCrossSection(G4int Z,
							     G4double incidentEnergy,
							     G4double mass,
							     G4double deltaEnergy,
							     G4bool testFlag) const
{

  deltaEnergy = 0;
  testFlag = 0;

  std::vector<G4double> crossSections;

  crossSections.push_back( paulShellK->CalculateKCrossSection(Z, mass, incidentEnergy) );
  
  G4Proton* aProtone = G4Proton::Proton();
  
  if (mass == aProtone->GetPDGMass() ) {
    
    crossSections.push_back( orlicShellLi->CalculateL1CrossSection(Z, incidentEnergy) );
    crossSections.push_back( orlicShellLi->CalculateL2CrossSection(Z, incidentEnergy) );
    crossSections.push_back( orlicShellLi->CalculateL3CrossSection(Z, incidentEnergy) );
  }
  
  return crossSections;

}




std::vector<G4double> G4empCrossSection::Probabilities(G4int Z,
							   G4double incidentEnergy,
							   G4double mass,
							   G4double deltaEnergy) const
{
  
std::vector<G4double> crossSections = GetCrossSection(Z, incidentEnergy, mass, deltaEnergy);

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



