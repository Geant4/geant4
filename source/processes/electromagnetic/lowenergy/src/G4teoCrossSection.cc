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
//$Id: G4teoCrossSection.cc,v 1.2 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Haifa Ben Abdelouahed
//         
//
// History:
// -----------
//  21 Apr 2009   ALF   1st implementation
//
// -------------------------------------------------------------------
// Class description:
// Low Energy Electromagnetic Physics, Cross section, p ionisation, K shell
// Further documentation available from http://www.ge.infn.it/geant4/lowE
// -------------------------------------------------------------------


#include "globals.hh"
#include "G4teoCrossSection.hh"
//#include "G4AtomicTransitionManager.hh"
//#include "G4NistManager.hh"
#include "G4Proton.hh"
//#include "G4Alpha.hh"
//#include <math.h>

G4teoCrossSection::G4teoCrossSection(G4String shellModel)
  :totalCS(0)
{ 

  if (shellModel == "ecpssr") {

  ecpssrShellK = new G4ecpssrKCrossSection();

  }
  orlicShellLi = new G4OrlicLiCrossSection();
}

G4teoCrossSection::~G4teoCrossSection()
{ 

  delete ecpssrShellK;
  delete orlicShellLi;

}

std::vector<G4double> G4teoCrossSection::GetCrossSection(G4int Z,
							     G4double incidentEnergy,
							     G4double mass,
							     G4double deltaEnergy,
							     G4bool testFlag) const
{

  std::vector<G4double> crossSections;

  crossSections.push_back( ecpssrShellK->CalculateCrossSection(Z, mass, incidentEnergy) );
  
  G4Proton* aProtone = G4Proton::Proton();
  
  if (mass == aProtone->GetPDGMass() ) {
    
    crossSections.push_back( orlicShellLi->CalculateL1CrossSection(Z, incidentEnergy) );
    crossSections.push_back( orlicShellLi->CalculateL2CrossSection(Z, incidentEnergy) );
    crossSections.push_back( orlicShellLi->CalculateL3CrossSection(Z, incidentEnergy) );
  }
  
  return crossSections;

}




std::vector<G4double> G4teoCrossSection::Probabilities(G4int Z,
							   G4double incidentEnergy,
							   G4double mass,
							   G4double deltaEnergy) const
{
  
std::vector<G4double> crossSections = GetCrossSection(Z, incidentEnergy, mass, deltaEnergy);

  for (G4int i=0; i<crossSections.size(); i++ ) {
    
    if (totalCS) {
      crossSections[i] = crossSections[i]/totalCS;
    }
    
  }

  return crossSections;

}


void G4teoCrossSection::SetTotalCS(G4double val){

  totalCS = val;

}



