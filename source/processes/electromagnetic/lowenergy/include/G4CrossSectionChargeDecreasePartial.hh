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
// $Id: G4CrossSectionChargeDecreasePartial.hh,v 1.1 2007/11/08 18:24:04 pia Exp $
// GEANT4 tag $Name: geant4-09-01 $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// Date         Name              Modification
// 28 Apr 2007  M.G. Pia          Created in compliance with design described in TNS paper
// 08 Nov 2007  MGP               Got code from S.I.; lowercase data member
//
// -------------------------------------------------------------------

// Class description:
// Geant4-DNA Cross total cross section for electron elastic scattering in water
// Reference: TNS Geant4-DNA paper
// S. Chauvie et al., Geant4 physics processes for microdosimetry simulation:
// design foundation and implementation of the first set of models,
// IEEE Trans. Nucl. Sci., vol. 54, no. 6, Dec. 2007.
// Reference for implementation model: NIM. 155, pp. 145-156, 1978
// Further documentation available from http://www.ge.infn.it/geant4/dna

// -------------------------------------------------------------------


#ifndef G4CROSSSECTIONCHARGEDECREASEPARTIAL_HH
#define G4CROSSSECTIONCHARGEDECREASEPARTIAL_HH 1
 
#include "globals.hh"

class G4ParticleDefinition;
 
class G4CrossSectionChargeDecreasePartial
{
public:
  
  G4CrossSectionChargeDecreasePartial();
  
  virtual ~G4CrossSectionChargeDecreasePartial();
  
  G4double CrossSection(G4double energy, G4int level, const G4ParticleDefinition* particle);

  G4double Sum(G4double energy, const G4ParticleDefinition* particle);

  G4int RandomSelect(G4double energy, const G4ParticleDefinition* particle);
  
  // Copy constructor and assignment operator to be added here
    
private:
   
   G4int numberOfPartialCrossSections[3]; // 3 is the particle type index

   G4double f0[2][3];
   G4double a0[2][3];
   G4double a1[2][3];
   G4double b0[2][3];
   G4double b1[2][3];
   G4double c0[2][3];
   G4double d0[2][3];
   G4double x0[2][3];
   G4double x1[2][3];

};

#endif
