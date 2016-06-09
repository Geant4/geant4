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
// Implementation of the HETC88 code into Geant4.
// Evaporation and De-excitation parts
// T. Lampen, Helsinki Institute of Physics, May-2000

#ifndef G4BEGammaDeexcitation_h
#define G4BEGammaDeexcitation_h 1

#include "globals.hh"

class G4BEGammaDeexcitation
{
public:
  G4BEGammaDeexcitation();
  virtual ~G4BEGammaDeexcitation();

  void setVerboseLevel( G4int verbose ); 

  void setNucleusA( G4int inputA );
  void setNucleusZ( G4int inputZ );
  void setExcitationEnergy( G4double inputE );

  G4DynamicParticle * emit();

private:  
  G4double sampleKineticEnergy();
  G4int verboseLevel;
  G4int nucleusA;
  G4int nucleusZ;
  G4double excitationEnergy;
  void isotropicCosines( G4double&,
			 G4double&,
			 G4double& );
};


#endif
