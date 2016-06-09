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
// $Id: G4CollisionPN.cc,v 1.2 2003/11/03 17:53:28 hpw Exp $ //
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//      For information related to this code contact:
//
//      File name:     G4CollisionNN
//
//      Author:        Maria Grazia Pia
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------


#include "globals.hh"
#include "G4CollisionPN.hh"
#include "G4CollisionComposite.hh"
#include "G4VCollision.hh"
#include "G4CollisionVector.hh"
#include "G4KineticTrack.hh"
#include "G4VCrossSectionSource.hh"
#include "G4XnpTotal.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4CollisionnpElastic.hh"
#include "G4CollisionNNToNDelta.hh"


G4CollisionPN::G4CollisionPN()
{ 

  crossSectionSource = new G4XnpTotal();

  // Components vector to be filled
  components = new G4CollisionVector;
  components->push_back(new G4CollisionnpElastic());
  components->push_back(new G4CollisionNNToNDelta());

  // Subtype of interacting particles
  G4String subType1 = G4Proton::ProtonDefinition()->GetParticleSubType();
  G4String subType2 = G4Neutron::NeutronDefinition()->GetParticleSubType();

  colliders1.push_back(subType1);
  colliders2.push_back(subType2);

}


G4CollisionPN::~G4CollisionPN()
{ 
  delete crossSectionSource;
  crossSectionSource = 0;

  if (components != 0)
    {
      G4int nComponents = components->size();
      for (G4int i = 0; i < nComponents; i++)
	{
	  G4CollisionPtr componentPtr = (*components)[i];
	  G4VCollision* component = componentPtr();
	  delete component;
	  component = 0;
	}
    }
  delete components;
  components = 0; 
}


const std::vector<G4String>& G4CollisionPN::GetListOfColliders(G4int whichOne) const
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
	  throw G4HadronicException(__FILE__, __LINE__, "G4CollisionNN::GetListOfColliders - Argument outside valid range"); 
	  return colliders1;
	}
    }
}


