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
// $Id: G4CrossSectionExcitationMillerGreen.hh,v 1.1 2007-05-02 17:18:48 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// Reference for implementation model: NIM. 155, pp. 145-156, 1978
// Further documentation available from http://www.ge.infn.it/geant4/dna

// -------------------------------------------------------------------


#ifndef G4ECROSSSECTIONEXCITATIONEMFIETZOGLOU_HH
#define G4ECROSSSECTIONEXCITATIONEMFIETZOGLOU_HH 1
 
#include "globals.hh"
#include "G4Track.hh"
 
 class G4CrossSectionExcitationMillerGreen
 {
  public:

   G4CrossSectionExcitationMillerGreen();
 
   virtual ~G4CrossSectionExcitationMillerGreen();
   
   G4double CrossSection(const G4Track&);

   // Copy constructor and assignment operator to be added here

 private:
   
   G4double EnergyConstant(G4int excitationLevel);

   G4double PartialCrossSection(G4double k, 
				G4int z, 
				G4int excitationLevel, 
				const G4ParticleDefinition* particle);

G4double S_1s(G4double t, 
	      G4double energyTransferred, 
	      G4double slaterEffectiveCharge, 
	      G4double shellNumber);

G4double S_2s(G4double t, 
	      G4double energyTransferred, 
	      G4double slaterEffectiveCharge, 
	      G4double shellNumber);


G4double S_2p(G4double t, 
	      G4double energyTransferred, 
	      G4double slaterEffectiveCharge, 
	      G4double shellNumber);


G4double R(G4double t, 
	   G4double energyTransferred, 
	   G4double slaterEffectiveCharge, 
	   G4double shellNumber);

   //  G4int RandomizePartialCrossSection(G4double k, G4int z);
 
   G4String name;  
   G4double lowEnergyLimit;
   G4double highEnergyLimit;

   G4double kineticEnergyCorrection;

   G4double slaterEffectiveCharge[3];
   G4double sCoefficient[3];
 };


#endif
