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
//         
//
// History:
// -----------
//  21 Apr 2009   ALF   1st implementation
//  29 Apr 2009   ALF Updated Desing for Integration
//  15 Mar 2011   ALF introduced the usage of G4AtomicShellEnumerator
//  20 Oct 2011   ALF updated to take into account ECPSSR Form Factor
//  09 Mar 2012   LP  update methods
//  09 Mar 2012   ALF  update for M-shells Simulation
//

#include "globals.hh"
#include "G4teoCrossSection.hh"
#include "G4Proton.hh"
#include "G4ecpssrBaseKxsModel.hh"
#include "G4ecpssrBaseLixsModel.hh"

#include "G4ecpssrFormFactorKxsModel.hh"
#include "G4ecpssrFormFactorLixsModel.hh"
#include "G4ecpssrFormFactorMixsModel.hh"

G4teoCrossSection::G4teoCrossSection(const G4String& nam)
  :G4VhShellCrossSection(nam),totalCS(0.0)
{ 
  ecpssrShellMi = nullptr;
  if (nam == "ECPSSR_Analytical") 
    {
      ecpssrShellK  = new G4ecpssrBaseKxsModel();  
      ecpssrShellLi = new G4ecpssrBaseLixsModel();      
    }
  else if (nam == "ECPSSR_FormFactor")
    {
      ecpssrShellK  = new G4ecpssrFormFactorKxsModel();  
      ecpssrShellLi = new G4ecpssrFormFactorLixsModel();
      ecpssrShellMi = new G4ecpssrFormFactorMixsModel();
    }
  else 
    { 
      G4cout << "G4teoCrossSection::G4teoCrossSection: ERROR " 
	     << " in cross section name ECPSSR_Analytical is used"
	     << G4endl;
      ecpssrShellK  = new G4ecpssrBaseKxsModel();  
      ecpssrShellLi = new G4ecpssrBaseLixsModel();      
    }
}

G4teoCrossSection::~G4teoCrossSection()
{ 
  delete ecpssrShellK;
  delete ecpssrShellLi;
  delete ecpssrShellMi;
}

std::vector<G4double> G4teoCrossSection::GetCrossSection(G4int Z,
							 G4double incidentEnergy,
							 G4double mass,
							 G4double,
							 const G4Material*)
{
  std::vector<G4double> crossSections;

  crossSections.push_back( ecpssrShellK->CalculateCrossSection(Z, mass, incidentEnergy) );
  
  crossSections.push_back( ecpssrShellLi->CalculateL1CrossSection(Z, mass, incidentEnergy) );
  crossSections.push_back( ecpssrShellLi->CalculateL2CrossSection(Z, mass, incidentEnergy) );
  crossSections.push_back( ecpssrShellLi->CalculateL3CrossSection(Z, mass, incidentEnergy) );

  if (ecpssrShellMi) {

    crossSections.push_back( ecpssrShellMi->CalculateM1CrossSection(Z, mass, incidentEnergy) );
    crossSections.push_back( ecpssrShellMi->CalculateM2CrossSection(Z, mass, incidentEnergy) );
    crossSections.push_back( ecpssrShellMi->CalculateM3CrossSection(Z, mass, incidentEnergy) );
    crossSections.push_back( ecpssrShellMi->CalculateM4CrossSection(Z, mass, incidentEnergy) );
    crossSections.push_back( ecpssrShellMi->CalculateM5CrossSection(Z, mass, incidentEnergy) );

  }


  return crossSections;
}

G4double G4teoCrossSection::CrossSection(G4int Z, G4AtomicShellEnumerator shell,
					 G4double incidentEnergy,
					 G4double mass,
					 const G4Material*)
{
  G4double res = 0.0;
  if(shell > 3 && !ecpssrShellMi) {
    return res; 
  } 

  else if(shell > 8) { 
    return res;
  } 
  
  else if(fKShell  == shell) 
    { 
      res = ecpssrShellK->CalculateCrossSection(Z, mass, incidentEnergy);
    } 
  
  else if(fL1Shell == shell) 
    { 
      res = ecpssrShellLi->CalculateL1CrossSection(Z, mass, incidentEnergy);
    } 
  
  else if(fL2Shell == shell) 
    { 
      res = ecpssrShellLi->CalculateL2CrossSection(Z, mass, incidentEnergy);
    } 
  
  else if(fL3Shell == shell) 
    { 
      res = ecpssrShellLi->CalculateL3CrossSection(Z, mass, incidentEnergy);
    }

  else if(fM1Shell == shell) 
    { 
      res = ecpssrShellMi->CalculateM1CrossSection(Z, mass, incidentEnergy);
    }

  else if(fM2Shell == shell) 
    { 
      res = ecpssrShellMi->CalculateM2CrossSection(Z, mass, incidentEnergy);
    }

  else if(fM3Shell == shell) 
    { 
      res = ecpssrShellMi->CalculateM3CrossSection(Z, mass, incidentEnergy);
    }

  else if(fM4Shell == shell) 
    { 
      res = ecpssrShellMi->CalculateM4CrossSection(Z, mass, incidentEnergy);
    }

  else if(fM5Shell == shell) 
    { 
      res = ecpssrShellMi->CalculateM5CrossSection(Z, mass, incidentEnergy);
    }
  return res;
}

std::vector<G4double> G4teoCrossSection::Probabilities(G4int Z,
						       G4double incidentEnergy,
						       G4double mass,
						       G4double deltaEnergy,
						       const G4Material*) 
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



