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
// $Id: G4CascadeInterface.hh,v 1.16 2010-05-21 18:07:30 mkelsey Exp $
// Defines an interface to Bertini (BERT) cascade
// based on INUCL  intra-nuclear transport.models 
// with bullet hadron energy ~< 10 GeV
//
// 20100405  M. Kelsey -- Fix constness of op== and op!=
// 20100519  M. Kelsey -- Remove Collider data members

#ifndef G4CASCADEINTERFACE_H
#define G4CASCADEINTERFACE_H 1

#include "G4VIntraNuclearTransportModel.hh"
#include "G4FragmentVector.hh"
#include "G4KineticTrackVector.hh"
#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4ParticleChange.hh"
#include "G4ReactionProduct.hh"
#include "G4ReactionProductVector.hh"


class G4CascadeInterface : public G4VIntraNuclearTransportModel {

public:
  G4CascadeInterface(const G4String& name = "Bertini Cascade");

  virtual ~G4CascadeInterface();

  G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries,
				     G4V3DNucleus* theNucleus);
  
  G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
				 G4Nucleus& theNucleus); 

private:
  G4int operator==(const G4CascadeInterface& right) const {
    return (this == &right);
  }

  G4int operator!=(const G4CascadeInterface& right) const {
    return (this != &right);
  }

  G4int verboseLevel;

private:
  G4HadFinalState theResult;  
};

#endif // G4CASCADEINTERFACE_H
