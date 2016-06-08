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
// $Id: G4PVHitsCollection.cc,v 1.9 2001/07/11 10:02:15 gunter Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

// G4PVHitsCollection

#include <assert.h>
#include "G4PVHitsCollection.hh"
#include "G4PersistentHitMan.hh"
#include "G4PHCofThisEvent.hh"

G4PVHitsCollection::G4PVHitsCollection(G4String detName,G4String colNam)
 : pcollectionName(colNam), pSDname(detName)
{;}

G4PVHitsCollection::~G4PVHitsCollection()
{;}

int G4PVHitsCollection::operator==(const G4PVHitsCollection &right) const
{ 
  return ((pcollectionName==right.pcollectionName)
        &&(pSDname==right.pSDname));
}

