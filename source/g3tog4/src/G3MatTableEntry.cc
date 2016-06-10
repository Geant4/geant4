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
// $Id: G3MatTableEntry.cc 67982 2013-03-13 10:36:03Z gcosmo $
//
// by I.Hrivnacova, 27 Sep 99

#include "G3MatTableEntry.hh"

#include "G4Material.hh"

G3MatTableEntry::G3MatTableEntry(G4int id, G4Material* material)
  : fID(id),
    fMaterial(material)
{}

G3MatTableEntry::G3MatTableEntry(const G3MatTableEntry& right)
  : fID(right.GetID()),
    fMaterial(right.GetMaterial())
{}    

G3MatTableEntry::~G3MatTableEntry()
{}

G3MatTableEntry& G3MatTableEntry::operator=(const G3MatTableEntry& right)
{
  if (&right == this)  { return *this; }
  fID = right.GetID();
  fMaterial = right.GetMaterial();     
  return *this;
}

G4int G3MatTableEntry::operator==(const G3MatTableEntry& right) const
{ 
  if (fID == right.GetID()) 
    return 1;
  else
    return 0;
}

G4int G3MatTableEntry::operator!=(const G3MatTableEntry& right) const
{ 
  if (*this == right) 
    return 0;
  else
    return 1;
}

