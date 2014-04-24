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
// $Id: G4CollisionNNElastic.cc,v 1.3 2006-06-29 20:37:32 gunter Exp $ //

#include "globals.hh"
#include "G4CollisionNNElastic.hh"
#include "G4KineticTrack.hh"
#include "G4VCrossSectionSource.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4XNNElastic.hh"
#include "G4AngularDistribution.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4KineticTrackVector.hh"
#include "G4AngularDistributionPP.hh"   //  nn and pp scattering

G4CollisionNNElastic::G4CollisionNNElastic()
{ 
  // Subtype of interacting particles
  G4String subType1 = G4Proton::ProtonDefinition()->GetParticleSubType();
  G4String subType2 = G4Neutron::NeutronDefinition()->GetParticleSubType();

  colliders1.push_back(subType1);
  colliders2.push_back(subType2);

  angularDistribution = new G4AngularDistributionPP();
  crossSectionSource = new G4XNNElastic();
}


G4CollisionNNElastic::~G4CollisionNNElastic()
{ 
  delete angularDistribution;
  angularDistribution = 0;
  delete crossSectionSource;
  crossSectionSource = 0;
}


G4bool G4CollisionNNElastic::IsInCharge(const G4KineticTrack& trk1, 
					const G4KineticTrack& trk2) const
{ 
  G4bool isInCharge = false;

  const G4ParticleDefinition* def1 = trk1.GetDefinition();
  const G4ParticleDefinition* def2 = trk2.GetDefinition();
  
  if ( (def1 == G4Proton::ProtonDefinition() && 
	def2 == G4Proton::ProtonDefinition() )
       ||
       (def1 == G4Neutron::NeutronDefinition() && 
	def2 == G4Neutron::NeutronDefinition() ) )
    {
      isInCharge = true;
    }
  return isInCharge;
}


const std::vector<G4String>& G4CollisionNNElastic::GetListOfColliders(G4int whichOne) const
{
  if (whichOne == 1) {
      return colliders1;
  }else if (whichOne == 2) {
	 return colliders2;
  }
  throw G4HadronicException(__FILE__, __LINE__, "G4CollisionNNElastic::GetListOfColliders - Argument outside valid range");
}
