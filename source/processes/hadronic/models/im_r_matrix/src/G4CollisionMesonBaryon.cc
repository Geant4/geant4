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

#include "globals.hh"
#include "G4CollisionMesonBaryon.hh"
#include "G4CollisionMesonBaryonElastic.hh"
#include "G4CollisionComposite.hh"
#include "G4VCollision.hh"
#include "G4CollisionVector.hh"
#include "G4KineticTrack.hh"
#include "G4VCrossSectionSource.hh"
#include "G4XNNTotal.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionZero.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"

#include "G4CollisionMesonBaryonToResonance.hh"
#include "G4CollisionMesonBaryonElastic.hh"

G4CollisionMesonBaryon::G4CollisionMesonBaryon()
{ 
  G4CollisionComposite::AddComponent(new G4CollisionMesonBaryonToResonance()); 
  G4CollisionComposite::AddComponent(new G4CollisionMesonBaryonElastic()); 
}

