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
// CLASS DESCRIPTION
// G4CascadeInterface defines an interface to HETC and INUCL 
// models of an medium energy (~ 0.5 - 5 GeV) intra-nuclear transport.
// If you have any questions, please contact 
// package writer aatos.heikkinen@cern.ch.
// --------------------------------------------------------------------

// This file is based on G4CascadeInterface.hh
// Modifications by: Pekka Kaitaniemi (kaitanie@cc.helsinki.fi)
// Helsinki Institute of Physics

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


class G4ElasticCascadeInterface : public G4VIntraNuclearTransportModel {

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
