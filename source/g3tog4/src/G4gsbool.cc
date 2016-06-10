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
// $Id: G4gsbool.cc 67982 2013-03-13 10:36:03Z gcosmo $
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
    G4String text = "G4gsbool: '" + volName + "' has no VolTableEntry";
    G4Exception("G4gsbool()", "G3toG40012", FatalException, text);
    return;
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
