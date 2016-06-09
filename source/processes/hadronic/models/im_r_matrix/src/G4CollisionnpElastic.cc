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

  G4ParticleDefinition* def1 = trk1.GetDefinition();
  G4ParticleDefinition* def2 = trk2.GetDefinition();
  
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
  if (whichOne == 1) 
    {
      return colliders1;
    }
  else 
    {
      if (whichOne == 2) 
	{ return colliders2; }
      else 
	{
	  throw G4HadronicException(__FILE__, __LINE__, "G4CollisionnpElastic::GetListOfColliders - Argument outside valid range"); 
	  return colliders1;
	}
    }
}
