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
//
// Implementation of the HETC88 code into Geant4.
// Evaporation and De-excitation parts
// T. Lampen, Helsinki Institute of Physics, May-2000

#ifndef G4VEvaporation_h
#define G4VEvaporation_h 1

#include "globals.hh"
#include "G4LayeredNucleus.hh"
#include "G4BertiniEvaporationChannel.hh"
#include "G4ParticleChange.hh"   


class G4BertiniEvaporation 
{

public:
  G4BertiniEvaporation();
  ~G4BertiniEvaporation(); 

  G4VParticleChange * BreakItUp( G4LayeredNucleus & nucleus);
  void setVerboseLevel( const G4int verbose );
  
private:  
  G4int verboseLevel;
  vector< G4BertiniEvaporationChannel * > channelVector;
  void G4BertiniEvaporation::fillParticleChange( vector< G4DynamicParticle * > secondaryParticleVector,
						 G4ParticleChange * theParticleChange );
  void G4BertiniEvaporation::splitBe8( const G4double E, 
				       const G4ThreeVector boost,
				       vector< G4DynamicParticle * > & secondaryParticleVector);
  void G4BertiniEvaporation::isotropicCosines( G4double & u, G4double & v,G4double & w );
};

#endif
