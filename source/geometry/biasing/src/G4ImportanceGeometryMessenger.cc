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
// $Id: G4ImportanceGeometryMessenger.cc,v 1.2 2002-07-12 10:40:44 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ImportanceGeometryMessenger.cc
//
// ----------------------------------------------------------------------

#include "G4ImportanceGeometryMessenger.hh"
#include "G4ImportanceGeometryConstructor.hh"
#include "G4UIcmdWithAString.hh"

G4ImportanceGeometryMessenger::
G4ImportanceGeometryMessenger(G4ImportanceGeometryConstructor &igeo):
  fImpGeoConst(igeo)
{
  fSolidTypeCmd = new G4UIcmdWithAString("/imp/worldvolume/solid",
					 this);
}

void G4ImportanceGeometryMessenger::
SetNewValue(G4UIcommand * command, G4String newValue){
  
  if (command==fSolidTypeCmd) {
    fImpGeoConst.SetSolidType(newValue);
  }
  
}
