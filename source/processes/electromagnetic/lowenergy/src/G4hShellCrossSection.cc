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
//
// ---------------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:  G4hShellCrossSection   
//
// Author:     S. Dussoni and A. Mantero (Alfonso.Mantero@ge.infn.it)
// 
// History:
// -----------
// 23 Oct 2001 A. Mantero   1st implementation
// 24 Oct 2001 MGP          Cleaned up
// 30 Oct 2001 V.Ivanchenko Include formula (53)
// 07 Oct 2002 V.Ivanchenko Fix in formula (53)
// 22 Apr 2004 S.Saliceti   Add testFlag variable and GetCrossSection method
// ----------------------------------------------------------------------------

#include "globals.hh"
#include "G4hShellCrossSection.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4ParticleDefinition.hh"

G4hShellCrossSection::G4hShellCrossSection()
{ }


G4hShellCrossSection::~G4hShellCrossSection()
{ }



std::vector<G4double> G4hShellCrossSection::GetCrossSection(G4int Z,
						G4double incidentEnergy,
						G4double mass,
						G4double deltaEnergy,
						G4bool) const 
{
  return CalculateCrossSections(Z,incidentEnergy,mass,deltaEnergy,false);
}  


std::vector<G4double> G4hShellCrossSection::CalculateCrossSections(G4int Z, 
								   G4double incidentEnergy, 
								   G4double hMass, 
								   G4double deltaEnergy,
								   G4bool testFlag) const
{
 // Cross-sections for proton ionization calculated as in
 // "M. Gryzinski, Two-Particle Collisions. I. General Relations for 
 // Collisions in the Laboratory system, Phys.Rev. 138 A305"
 // Other reference papers are Gryzinski's "Paper I" and "Paper II"
 // V.Ivanchenko add only implementation of the formula (53) 
 // last factor neglected because it is 1 with a good accuracy

 const G4AtomicTransitionManager*  transitionManager = 
                             G4AtomicTransitionManager::Instance();

 size_t nShells = transitionManager->NumberOfShells(Z);

 // Vector that stores the calculated cross-sections for each shell:
 std::vector<G4double> crossSections;

 // In this loop we calculate cross-section for every shell in the atom
 for (size_t k=0;  k<nShells;  k++)
   {
     G4double bindingEnergy = transitionManager->Shell(Z,k)->BindingEnergy();
     G4double xEnergy = 0.5*bindingEnergy;
   G4double y = incidentEnergy*electron_mass_c2/(xEnergy*hMass);
     G4double dele = deltaEnergy + xEnergy;
     G4double x = electron_mass_c2/dele;

     G4double aCrossSection = (dele/(xEnergy*(1. + 1./y)) 
                            + 4.*std::log(2.7 + std::sqrt(y))/3.) * x*x*x;    

      // Fill the vector of cross sections with the value just calculated
     crossSections.push_back(aCrossSection);

      if (testFlag) 
 	{
	  G4cout <<"Element: " <<Z<<" Shell: "<<k<<" Particle Energy(MeV): "<<incidentEnergy/MeV<<G4endl;
         G4cout <<"Delta Ray Energy(MeV): "<< deltaEnergy/MeV<<" Binding Energy(MeV): "<<bindingEnergy/MeV<<G4endl;
	  G4cout <<"Cross Section: "<<aCrossSection/barn<<" barns"<< G4endl;
	}
    }

 return crossSections;
}



// Normalization of relative cross-sections to 1

std::vector<G4double> G4hShellCrossSection::Probabilities(
                                              G4int Z, 
			                      G4double incidentEnergy, 
					      G4double hMass, 
					      G4double deltaEnergy
					      ) const
{  
  std::vector<G4double> crossSection = CalculateCrossSections(Z, incidentEnergy, hMass, deltaEnergy,true);
  G4double nShells = crossSection.size();

  // Partial and total cross-section used for normalization of crossSections:
  G4double totalCrossSection = 0.;
  
  // Calculation of total cross-section
  for (size_t j=0; j<nShells; j++)
    {      
      totalCrossSection += crossSection[j];  // equivalent to totalCrossSection = totalCrossSection + aCrossSection  
    }

  // Calculation of normal cross section
  for(size_t i=0; i<nShells; i++)     
    {
      crossSection[i] = crossSection[i] / totalCrossSection;
    }
  // Returns the normalized vector
  return crossSection;
}
  
  

