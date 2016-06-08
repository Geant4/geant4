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
// $Id: G4PVDigitsCollection.cc,v 1.4 2001/07/11 10:02:13 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//

// G4PVDigitsCollection

#include <assert.h>
#include "G4PVDigitsCollection.hh"
#include "G4PersistentDigitMan.hh"
#include "G4PDCofThisEvent.hh"

G4PVDigitsCollection::G4PVDigitsCollection(G4String detName,G4String colNam)
 : pcollectionName(colNam), pDMname(detName)
{;}

G4PVDigitsCollection::~G4PVDigitsCollection()
{;}

int G4PVDigitsCollection::operator==(const G4PVDigitsCollection &right) const
{ 
  return ((pcollectionName==right.pcollectionName)
        &&(pDMname==right.pDMname));
}

