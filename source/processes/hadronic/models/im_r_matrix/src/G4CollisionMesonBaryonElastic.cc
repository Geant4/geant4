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
