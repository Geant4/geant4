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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:  G4hShellCrossSectionExp.cc   
//
// Author:     Simona Saliceti (simona.saliceti@ge.infn.it)
// 
// History:
// -----------
// From 23 Oct 2001 A. Mantero G4hShellCrossSection
// 30/03/2004 Simona Saliceti 1st implementation
// -------------------------------------------------------------------
// Class Description: 
// Empiric Model for shell cross sections in proton ionisation
// -------------------------------------------------------------------
// $Id: G4hShellCrossSectionExp.cc,v 1.7 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "globals.hh"
#include <vector>
#include "G4hShellCrossSectionExp.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4Electron.hh"
#include "G4hShellCrossSectionExpData.hh"
#include "G4Proton.hh"
#include "G4ParticleDefinition.hh"

G4hShellCrossSectionExp::G4hShellCrossSectionExp()
{
  kShellData = new G4hShellCrossSectionExpData();
}

G4hShellCrossSectionExp::~G4hShellCrossSectionExp()
{ }

std::vector<G4double> G4hShellCrossSectionExp::GetCrossSection(G4int Z,
							       G4double incidentEnergy,
							       G4double mass,
							       G4double deltaEnergy,
							       G4bool testFlag) const
{
  mass = 0.0;
  deltaEnergy = 0.0;

  std::vector<G4double> aCrossSection; 
 
  // Fill the vector of cross sections with the value just calculated
  aCrossSection.push_back(GetCrossSectionExp(Z,incidentEnergy));
 
  if (testFlag) 
    {
      G4cout <<"Element: " <<Z<<" Particle Energy: "<<incidentEnergy/MeV<<" MeV" <<G4endl;
      G4cout <<"Cross Section: "<<aCrossSection[0]/barn<<" barns"<< G4endl;
    }
  return aCrossSection;
}

//This function calculated the cross section with the Empiric model
G4double G4hShellCrossSectionExp::GetCrossSectionExp(G4int Z, 
						     G4double incidentEnergy) const
{
  // Vector that stores the calculated cross-sections for each shell:
  G4double  crossSectionsInBarn = 0.0;
  G4double  crossSections = 0.0;

  std::vector<G4double>* parVec = kShellData->GetParam(Z);
  std::vector<G4double>::iterator i = (*parVec).begin();
       
  G4double a = *i;
  G4double b = *(i+1);
  G4double c = *(i+2);

  G4double incidentEnergyInMeV = incidentEnergy/MeV;

  if(Z<26 && Z>=6 && Z!=17)
    {
      crossSectionsInBarn = a*(std::pow(b,(1./incidentEnergyInMeV)))*(std::pow(incidentEnergyInMeV,c));
    }          
  else if((Z<66 && Z>=26) || Z==17)
    {
      crossSectionsInBarn = std::exp(a+(b/incidentEnergyInMeV)+(c*std::log(incidentEnergyInMeV)));
    }
  else if(Z<=92 && Z>=66)
    {	 
      crossSectionsInBarn = (std::pow(incidentEnergyInMeV,a))*std::exp(b-(c*incidentEnergyInMeV));  
    }
  else
    {
      G4cout << "Error: there is not that Z" << G4endl;
    }     

  crossSections = crossSectionsInBarn*barn;
  return crossSections;
}

// This function gives the atomic cross section of k shell only
void G4hShellCrossSectionExp::SetTotalCS(G4double value)
{
  atomTotalCrossSection = value;
}

//A new implementation of Probability to calculate the cross section probability for k shell only
 std::vector<G4double> G4hShellCrossSectionExp::Probabilities(
							      G4int Z, 
							      G4double incidentEnergy, 
							      G4double hMass, 
							      G4double deltaEnergy
							      ) const
{  
  hMass = 0.0;
  deltaEnergy = 0.0;
  
  std::vector<G4double> kProbability;
  kProbability.push_back(GetCrossSectionExp(Z,incidentEnergy)/atomTotalCrossSection);
  // ---- MGP ---- Next line corrected to kProbability[0] instead of [1], which is not initialized!
  kProbability.push_back(1 - kProbability[0]);

  return kProbability;
}
  
  

