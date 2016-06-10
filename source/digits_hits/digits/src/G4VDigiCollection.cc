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
// $Id: G4VDigiCollection.cc 67992 2013-03-13 10:59:57Z gcosmo $
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

