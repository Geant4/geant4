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
// $Id: G4FinalStateIonisationRudd.hh,v 1.2 2007/11/09 16:30:56 pia Exp $
// GEANT4 tag $Name: geant4-09-01 $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// Date         Name              Modification
// 28 Apr 2007  M.G. Pia          Created in compliance with design described in TNS paper
//
// -------------------------------------------------------------------

// Class description:
// Geant4-DNA Cross total cross section for electron elastic scattering in water
// Reference: TNS Geant4-DNA paper
// Reference: TNS Geant4-DNA paper
// S. Chauvie et al., Geant4 physics processes for microdosimetry simulation:
// design foundation and implementation of the first set of models,
// IEEE Trans. Nucl. Sci., vol. 54, no. 6, Dec. 2007.
// Reference for implementation model: NIM. 155, pp. 145-156, 1978
// Further documentation available from http://www.ge.infn.it/geant4/dna

// -------------------------------------------------------------------


#ifndef G4FINALSTATEIONISATIONRUDD_HH
#define G4FINALSTATEIONISATIONRUDD_HH 1
 
#include "globals.hh"
#include "G4FinalStateProduct.hh"
#include "G4WaterIonisationStructure.hh"
#include "G4CrossSectionIonisationRuddPartial.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;

class G4FinalStateIonisationRudd
{
public:
   
  G4FinalStateIonisationRudd();
   
  ~G4FinalStateIonisationRudd();
   
  const G4FinalStateProduct& GenerateFinalState(const G4Track& track, const G4Step& step);
   
private:
   
  // Copy constructor and assignment operator to be added here
   
  G4String name;  
  G4double lowEnergyLimitDefault;
  G4double highEnergyLimitDefault;
  std::map<G4String,G4double,std::less<G4String> > lowEnergyLimit;
  std::map<G4String,G4double,std::less<G4String> > highEnergyLimit;

  G4FinalStateProduct product;

  G4WaterIonisationStructure waterStructure;

  G4CrossSectionIonisationRuddPartial cross;
   
  G4double RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition, 
					  G4double incomingParticleEnergy, 
					  G4int shell);

  void RandomizeEjectedElectronDirection(G4ParticleDefinition* particleDefinition, 
					 G4double incomingParticleEnergy, 
					 G4double outgoingParticleEnergy, 
					 G4double cosTheta, 
					 G4double phi);
   
  G4double  DifferentialCrossSection(G4ParticleDefinition* particleDefinition, 
				   G4double k, 
				   G4double energyTransfer, 
				   G4int shell);

  G4double CorrectionFactor(G4ParticleDefinition* particleDefinition, G4double k);

  G4double S_1s(G4double t, 
		G4double energyTransferred, 
		G4double slaterEffectiveChg, 
		G4double shellNumber);

  G4double S_2s(G4double t, 
		G4double energyTransferred, 
		G4double slaterEffectiveChg, 
		G4double shellNumber);


  G4double S_2p(G4double t, 
		G4double energyTransferred, 
		G4double slaterEffectiveChg, 
		G4double shellNumber);

  G4double R(G4double t, 
	     G4double energyTransferred, 
	     G4double slaterEffectiveChg, 
	     G4double shellNumber) ;

  G4double slaterEffectiveCharge[3];
  G4double sCoefficient[3];
};


#endif
