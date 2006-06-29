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

#ifndef G4GeneratorPrecompoundInterface_h
#define G4GeneratorPrecompoundInterface_h 1

#include "G4Fancy3DNucleus.hh"
#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4VIntraNuclearTransportModel.hh"
#include "G4KineticTrackVector.hh"
#include "G4FragmentVector.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"

// Class Description
// Trivial implementation of an intra-nuclear transport. It pworvides coupling
// of high energy generators with pre equilibrium decay models. 
// To be used in your physics list in case you need this physics.
// Class Description - End

class G4GeneratorPrecompoundInterface : public G4VIntraNuclearTransportModel 
{
public:
   G4GeneratorPrecompoundInterface(){}      
   ~G4GeneratorPrecompoundInterface(){}

private:
   G4int operator==(G4GeneratorPrecompoundInterface& right) {return (this == &right);}
   G4int operator!=(G4GeneratorPrecompoundInterface& right) {return (this != &right);}
      
public:
   G4HadFinalState * ApplyYourself(const G4HadProjectile &aTrack, G4Nucleus &targetNucleus );
   G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus);


private:   
};

#endif // G4GeneratorPrecompoundInterface_h


