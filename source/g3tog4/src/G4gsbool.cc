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
// $Id: G4gsbool.cc,v 1.1 2001-11-08 16:08:01 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, 13.10.01

#include "G3G4Interface.hh"
#include "G3VolTable.hh"
#include "globals.hh"

void G4gsbool(G4String volName, G4String manyVolName)
{
  // find VTEs
  G3VolTableEntry* vte = G3Vol.GetVTE(volName);
  G3VolTableEntry* manyVTE = G3Vol.GetVTE(manyVolName);

  if (vte == 0) {
    G4Exception("G4gsbool: '" + volName + "' has no VolTableEntry");
  } 
  else if (manyVTE == 0) {
    // warning
    G4cerr << "G4gsbool: '" << manyVolName << "' has no VolTableEntry." 
           << G4endl
           << "          Specified overlap will be ignored."
	   << G4endl;
    return;	   
  } 
  else { 
    manyVTE->AddOverlap(vte);
  }   
}
