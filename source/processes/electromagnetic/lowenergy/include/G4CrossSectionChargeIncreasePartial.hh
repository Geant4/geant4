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
// $Id: G4CrossSectionChargeIncreasePartial.hh,v 1.2 2008-07-14 20:47:33 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4CROSSSECTIONCHARGEIncreasePARTIAL_HH
#define G4CROSSSECTIONCHARGEIncreasePARTIAL_HH 1
 
#include "G4DNAGenericIonsManager.hh"
#include "Randomize.hh"
 
class G4CrossSectionChargeIncreasePartial
{
public:
  
  G4CrossSectionChargeIncreasePartial();
  
  virtual ~G4CrossSectionChargeIncreasePartial();
  
  G4double CrossSection(G4double energy, G4int level, const G4ParticleDefinition* particle);

  G4double Sum(G4double energy, const G4ParticleDefinition* particle);

  G4int RandomSelect(G4double energy, const G4ParticleDefinition* particle);
  
private:
   
  G4int numberOfPartialCrossSections[2]; // 2 is the particle type index

  G4double f0[2][2];
  G4double a0[2][2];
  G4double a1[2][2];
  G4double b0[2][2];
  G4double b1[2][2];
  G4double c0[2][2];
  G4double d0[2][2];
  G4double x0[2][2];
  G4double x1[2][2];

};

#endif
