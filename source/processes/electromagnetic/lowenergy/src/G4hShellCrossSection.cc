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
// -------------------------------------------------------------------
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
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4hShellCrossSection.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4Electron.hh"

G4hShellCrossSection::G4hShellCrossSection()
{ }


G4hShellCrossSection::~G4hShellCrossSection()
{ }


G4std::vector<G4double> G4hShellCrossSection::Probabilities(G4int Z, 
			                                    G4double eKin, 
							    G4double mass, 
							    G4double p) const
{
  // Cross-sections for proton ionization calculated as in
  // "M. Gryzinski, Two-Particle Collisions. I. General Relations for 
  // Collisions in the Laboratory system, Phys.Rev. 138 A305"
  // Other reference papers are Gryzinski's "Paper I" and "Paper II"

  G4AtomicTransitionManager*  transitionManager = 
                              G4AtomicTransitionManager::Instance();

  size_t nShells = transitionManager->NumberOfShells(Z);

  //  G4cout << "Z= " << Z << " nShells= " << nShells << G4endl;

  // Vector that stores the calculated cross-sections for each shell:
  G4std::vector<G4double> crossSections;

  // Partial and total cross-section used for normalization of crossSections:
  G4double aCrossSection = 0.;
  G4double totalCrossSection = 0.;

  // k2 in Gryzinski's "Paper II"
  G4double incFactor = eKin/mass;	
  // Relativistically corrected velocity
  G4double beta = p * c_light / (sqrt(mass * mass * c_squared + p * p));
  // gamma
  G4double gamma =1. / sqrt(1. - beta*beta / c_squared);
  G4double elMass = G4Electron::ElectronDefinition()->GetPDGMass();
  G4double inverseMaxETransfer = eKin * 
    ( 1. + 2. * gamma*elMass / mass + ((elMass/mass)*(elMass/mass)) /
     (2. * elMass * (beta * beta / (c_squared - beta*beta) ) ) );

  // In this loop we calculate cross-section for every shell in the atom
  for (size_t k=0;  k<nShells;  k++)
    {
      G4double bindingEnergy = transitionManager->Shell(Z,k)->BindingEnergy();
    
    // Electron kinematics derived from its shell
      G4double elVelocity = sqrt(c_squared*(1. - elMass*elMass / (bindingEnergy*bindingEnergy)));	

      // k1 in Gryzinski's "Paper II"
      G4double elFactor = bindingEnergy / elMass;	

      // Begin of calculation of shell cross-section
      G4double eFactor1 = 1.+ elFactor;
      G4double eFactor2 = 2. + elFactor;
      G4double iFactor1 = 1.+ incFactor;
      G4double iFactor2 = 2. + incFactor;
      // - MGP - Probably there is something wrong in the calculation of velocityFunction
      // - MGP - To be checked in the original paper
      G4double velocityFunction = elFactor / incFactor *
	eFactor2 / iFactor2 * (iFactor1 * iFactor1) * eFactor1 / (eFactor1 * eFactor1) + 
	elFactor / incFactor * 	pow(eFactor2/iFactor2, 1.5);
      
      G4double eV2 = elVelocity * elVelocity;
      G4double beta2 = beta * beta;
      G4double secondFunction = beta2 / (beta2 + eV2) +
	2./3. * (1.+inverseMaxETransfer) * log(2.7+beta/elVelocity) * 
	(1.-inverseMaxETransfer) *
	pow((1.-inverseMaxETransfer),(1. + beta2 /eV2));
    
    aCrossSection = velocityFunction * secondFunction / (bindingEnergy * bindingEnergy);
    
    // VI only for test
    //  G4cout << "cross[" << k << "]= " << aCrossSection << G4endl;
    aCrossSection = 1.;
    // Calculation of total cross-section
    totalCrossSection += aCrossSection;
    
    // Fill the vector of cross sections with the value just calculated
    crossSections.push_back(aCrossSection);    
    }

  // Normalization of relative cross-sections to 1
  for (size_t j=0; j<nShells; j++)
    {
      crossSections[j] = crossSections[j] / totalCrossSection;
    }
  
  // Returns the normalized vector 
  return crossSections;
}
