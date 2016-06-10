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
//hp

#include "globals.hh"
#include "G4CollisionnpElastic.hh"
#include "G4KineticTrack.hh"
#include "G4VCrossSectionSource.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4XnpElastic.hh"
#include "G4AngularDistribution.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4KineticTrackVector.hh"
#include "G4AngularDistributionNP.hh"   //  np scattering
// hpw ok.
G4CollisionnpElastic::G4CollisionnpElastic()
{ 
  // Subtype of interacting particles
  G4String subType1 = G4Proton::ProtonDefinition()->GetParticleSubType();
  G4String subType2 = G4Neutron::NeutronDefinition()->GetParticleSubType();

  colliders1.push_back(subType1);
  colliders2.push_back(subType2);

//  angularDistribution = new G4AngularDistribution(true);
  angularDistribution = new G4AngularDistributionNP();
  crossSectionSource = new G4XnpElastic();
}


G4CollisionnpElastic::~G4CollisionnpElastic()
{ 
  delete angularDistribution;
  delete crossSectionSource;
}


G4bool G4CollisionnpElastic::IsInCharge(const G4KineticTrack& trk1, 
					const G4KineticTrack& trk2) const
{ 
  G4bool isInCharge = false;

  const G4ParticleDefinition* def1 = trk1.GetDefinition();
  const G4ParticleDefinition* def2 = trk2.GetDefinition();
  
  if ( (def1 == G4Neutron::NeutronDefinition() && 
	def2 == G4Proton::ProtonDefinition() )
       ||
       (def1 == G4Proton::ProtonDefinition() && 
	def2 == G4Neutron::NeutronDefinition() ) )
    {
      isInCharge = true;
    }
  return isInCharge;
}


const std::vector<G4String>& G4CollisionnpElastic::GetListOfColliders(G4int whichOne) const
{
	if (whichOne == 1)  {
		return colliders1;
	} else if (whichOne == 2) {
		return colliders2;
	}

	throw G4HadronicException(__FILE__, __LINE__, "G4CollisionnpElastic::GetListOfColliders - Argument outside valid range");
}
