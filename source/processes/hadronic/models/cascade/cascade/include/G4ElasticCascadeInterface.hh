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
// CLASS DESCRIPTION
// G4ElasticCascadeInterface defines an interface to INUCL 
// models of an medium energy (~ 0.5 - 10 GeV) intra-nuclear transport.
// Elastic reaction is forced in this interface
// If you have any questions, please contact 
// package writer aatos.heikkinen@cern.ch., 
// Also coded by Pekka Kaitataniemi, Helsinki Institute of Physics
// --------------------------------------------------------------------
#ifndef G4ELASTICCASCADEINTERFACE_H
#define G4ELASTICCASCADEINTERFACE_H 1

#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4VIntraNuclearTransportModel.hh"
#include "G4KineticTrackVector.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleChange.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"

#include "G4HadronicInteraction.hh"

//class G4CascadeInterface : public G4VIntraNuclearTransportModel {

//class G4ElasticCascadeInterface : public G4HadronElastic {

class G4ElasticCascadeInterface : public G4HadronicInteraction {
public:

  G4ElasticCascadeInterface();

  ~G4ElasticCascadeInterface(){
  }

  G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus);

  G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack, 
                                   G4Nucleus& theNucleus); 

private:

  G4int operator==(G4ElasticCascadeInterface& right) {
    return (this == &right);
  }

  G4int operator!=(G4ElasticCascadeInterface& right) {
    return (this != &right);
  }

  G4int verboseLevel;
private:
  G4HadFinalState theResult;  
  
};

#endif  //G4ELASTICCASCADEINTERFACE_H
