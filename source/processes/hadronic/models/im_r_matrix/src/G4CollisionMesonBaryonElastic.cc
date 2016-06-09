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
#include "G4CollisionMesonBaryonElastic.hh"
#include "G4AngularDistribution.hh"
#include "G4XMesonBaryonElastic.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionZero.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include <vector>
G4CollisionMesonBaryonElastic::G4CollisionMesonBaryonElastic()
{
  angularDistribution = new G4AngularDistribution(false);
  crossSectionSource = new G4XMesonBaryonElastic();
}

G4CollisionMesonBaryonElastic::~G4CollisionMesonBaryonElastic()
{ 
  delete angularDistribution;
  delete crossSectionSource;
}

G4bool G4CollisionMesonBaryonElastic::
 IsInCharge(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const
 {
   G4bool result = false;
   G4ParticleDefinition * p1 = trk1.GetDefinition();
   G4ParticleDefinition * p2 = trk2.GetDefinition();
   if(   (GetNumberOfPartons(p1) != 2 || GetNumberOfPartons(p2) != 3)
       ||(GetNumberOfPartons(p1) != 3 || GetNumberOfPartons(p2) != 2) ) 
   {
     result = false;
   }
   return result;
 }

G4String G4CollisionMesonBaryonElastic::
 GetName() const 
 { 
   return "Meson Baryon Elastic Collision"; 
 }

const std::vector<G4String>& G4CollisionMesonBaryonElastic::
 GetListOfColliders(G4int ) const
 {
   throw G4HadronicException(__FILE__, __LINE__, "Called G4CollisionMesonBaryonElastic::GetListOfColliders");
   return dummy;
 }  
