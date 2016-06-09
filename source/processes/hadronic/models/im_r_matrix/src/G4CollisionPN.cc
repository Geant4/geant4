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
#include "G4CollisionPN.hh"
#include "G4XnpTotal.hh"
#include "G4CollisionnpElastic.hh"
#include "G4CollisionNNToNDelta.hh"
#include "G4Pair.hh"

// J.P. Wellisch, Dec 2004.

typedef GROUP2(G4CollisionnpElastic, G4CollisionNNToNDelta) theChannels;

G4CollisionPN::G4CollisionPN()
{ 

  crossSectionSource = new G4XnpTotal();
  Register aR;
  G4ForEach<theChannels>::Apply(&aR, this);
}


G4CollisionPN::~G4CollisionPN()
{ 
  delete crossSectionSource;
  crossSectionSource = 0;
}


const std::vector<G4String>& G4CollisionPN::GetListOfColliders(G4int ) const
{
  throw G4HadronicException(__FILE__, __LINE__, "G4CollisionNN:: GetListOfColliders called");
  return colliders1;
}

