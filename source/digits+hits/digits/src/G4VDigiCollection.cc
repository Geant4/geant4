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
// $Id: G4VDigiCollection.cc,v 1.6 2001/07/13 15:00:15 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

// G4VDigiCollection

#include "G4VDigiCollection.hh"

G4VDigiCollection::G4VDigiCollection()
{
  collectionName = "Unknown";
  DMname = "Unknown";
}

G4VDigiCollection::G4VDigiCollection(G4String DMnam,G4String colNam)
{
  collectionName = colNam;
  DMname = DMnam;
}

G4VDigiCollection::~G4VDigiCollection()
{ ; }

G4int G4VDigiCollection::operator==(const G4VDigiCollection &right) const
{ 
  return ((collectionName==right.collectionName)
        &&(DMname==right.DMname));
}

void G4VDigiCollection::DrawAllDigi() 
{;}

void G4VDigiCollection::PrintAllDigi() 
{;}

