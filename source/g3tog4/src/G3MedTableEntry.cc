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
// $Id: G3MedTableEntry.cc 67982 2013-03-13 10:36:03Z gcosmo $
//
// by I.Hrivnacova, 27 Sep 99

#include "G3MedTableEntry.hh"

#include "G4Material.hh"
#include "G4MagneticField.hh"
#include "G4UserLimits.hh"

G3MedTableEntry::G3MedTableEntry(G4int id, G4Material* material, 
       G4MagneticField* field, G4UserLimits* limits, G4int isvol)
  : fID(id),
    fMaterial(material),
    fField(field),
    fLimits(limits),
    fISVOL(isvol)
{}

G3MedTableEntry::G3MedTableEntry(const G3MedTableEntry& right)
  : fID(right.GetID()),
    fMaterial(right.GetMaterial()),
    fField(right.GetField()),
    fLimits(right.GetLimits()),
    fISVOL(right.GetISVOL())    
{}    

G3MedTableEntry::~G3MedTableEntry()
{}

G3MedTableEntry& G3MedTableEntry::operator=(const G3MedTableEntry& right)
{ 
  if (&right == this)  { return *this; }
  fID = right.GetID();
  fMaterial = right.GetMaterial();     
  fField = right.GetField();
  fLimits = right.GetLimits();
  fISVOL = right.GetISVOL();  
  return *this;
}

G4int G3MedTableEntry::operator==(const G3MedTableEntry& right) const
{ 
  if (fID == right.GetID()) 
    return 1;
  else
    return 0;
}

G4int G3MedTableEntry::operator!=(const G3MedTableEntry& right) const
{ 
  if (*this == right) 
    return 0;
  else
    return 1;
}

