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
  throw G4HadronicException(__FILE__, __LINE__, "G4CollisionPN:: GetListOfColliders called");
  return colliders1;
}

