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
// $Id: G4VIsotopeTable.cc 66241 2012-12-13 18:34:42Z gunter $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
// **********************************************************************
//      New design using G4VIsotopeTable          30 Apr.. 2013 H.Kurashige

#include "G4VIsotopeTable.hh"

// ######################################################################
// ###                           IsotopeTable                         ###
// ######################################################################

#include "G4IsotopeProperty.hh"
#include "G4VIsotopeTable.hh"

G4VIsotopeTable::G4VIsotopeTable()
  : fName(""), verboseLevel(0)
{
}

G4VIsotopeTable::G4VIsotopeTable(const G4String& name)
  : fName(name), verboseLevel(0)
{
}

G4VIsotopeTable::G4VIsotopeTable(const G4VIsotopeTable & right)
  : fName(right.fName), verboseLevel(right.verboseLevel)
{
}

G4VIsotopeTable& G4VIsotopeTable::operator=(const G4VIsotopeTable & right)
{
  if (this != &right){
    fName = right.fName;
    verboseLevel = right.verboseLevel;
  }
  return *this;
}

G4VIsotopeTable::~G4VIsotopeTable()
{
}

G4IsotopeProperty* G4VIsotopeTable::GetIsotopeByIsoLvl(G4int Z, G4int A, G4int level)
{
  // temporal implementation
  if (level==0)  return GetIsotope(Z, A, 0.0);
  else           return 0;
}


void G4VIsotopeTable::DumpTable(G4int Zmin, G4int Zmax) 
{
  G4int Z, A;
  G4int lvl;
  const G4int MAX_LVL=9;
  for ( Z =Zmin; Z<=Zmax; Z++){
    for ( A= Z; A<=3*Z; A++){
      for ( lvl=0; lvl<=MAX_LVL; lvl++){
	G4IsotopeProperty* ptr = GetIsotope(Z,A,lvl);
	if (ptr!=0) ptr->DumpInfo();
      }
    }      
  }
}








