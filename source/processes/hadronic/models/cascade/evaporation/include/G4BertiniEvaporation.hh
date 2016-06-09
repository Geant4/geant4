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

#ifndef G4BertiniEvaporation_h
#define G4BertiniEvaporation_h 1

#include "globals.hh"
#include "G4LayeredNucleus.hh"
#include "G4BertiniEvaporationChannel.hh"
#include "G4ParticleChange.hh"   
#include "G4VEvaporation.hh"


class G4BertiniEvaporation : public G4VEvaporation
{

public:
  G4BertiniEvaporation();
  ~G4BertiniEvaporation(); 

  virtual G4FragmentVector * BreakItUp(const G4Fragment &theNucleus)
  {
    G4LayeredNucleus aNuc( theNucleus.GetA(), theNucleus.GetZ() );
    aNuc.AddExcitationEnergy(theNucleus.GetExcitationEnergy());
    return BreakItUp(aNuc);
  }
  G4FragmentVector * BreakItUp( G4LayeredNucleus & nucleus);
  void setVerboseLevel( const G4int verbose );
  
private:  
  G4int verboseLevel;
  std::vector< G4BertiniEvaporationChannel * > channelVector;
  void fillResult( std::vector< G4DynamicParticle * > secondaryParticleVector,
		   G4FragmentVector * aResult );
  void splitBe8( const G4double E, 
		 const G4ThreeVector boost,
		 std::vector< G4DynamicParticle * > & secondaryParticleVector);
  void isotropicCosines( G4double & u, G4double & v,G4double & w );
};

#endif
