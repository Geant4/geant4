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
// $Id: G4TDigiCollection.cc,v 1.4 2001/07/13 15:00:14 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-01 $
//

#include "G4TDigiCollection.hh"

G4Allocator<G4DigiCollection> aDCAllocator;

G4DigiCollection::G4DigiCollection()
{;}

G4DigiCollection::G4DigiCollection(G4String detName,G4String colNam)
: G4VDigiCollection(detName,colNam)
{;}

G4DigiCollection::~G4DigiCollection()
{;}

G4int G4DigiCollection::operator==(const G4DigiCollection &right) const
{ return (collectionName==right.collectionName); }

