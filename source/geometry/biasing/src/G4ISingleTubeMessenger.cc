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
// $Id: G4ISingleTubeMessenger.cc,v 1.2 2002-07-12 10:40:44 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ISingleTubeMessenger.cc
//
// ----------------------------------------------------------------------

#include "G4ISingleTubeMessenger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ITubeFactory.hh"

G4ISingleTubeMessenger::G4ISingleTubeMessenger(const G4String &cellname,
					       G4ITubeFactory *tfact) {

  fCellName = cellname;
  fITubeFactory = tfact;
  fZminIsSet = false;
  fZmaxIsSet = false;
  fIbasesIsSet = false;
  fIexpoIsSet = false;
  fICellCreated = false;
  
  G4String path("/imp/cell/" + fCellName + "/");

  fZminCmd = new G4UIcmdWithADoubleAndUnit(G4String(path + "zmin"),this);
  fZmaxCmd = new G4UIcmdWithADoubleAndUnit(G4String(path + "zmax"),this);
  fIbaseCmd = new G4UIcmdWithADouble(G4String(path + "ibase"),this);
  fIexpoCmd = new G4UIcmdWithADouble(G4String(path + "iexpo"),this);

}

void G4ISingleTubeMessenger::
SetNewValue(G4UIcommand * command, G4String newValue){

  if (command==fZminCmd) {
    fZmin = fZminCmd->GetNewDoubleValue(newValue);
    fZminIsSet = true;
  }
  if (command==fZmaxCmd) {
    fZmax = fZmaxCmd->GetNewDoubleValue(newValue);
    fZmaxIsSet = true;
  }
  if (command==fIbaseCmd) {
    fIbase = fIbaseCmd->GetNewDoubleValue(newValue);
    fIbasesIsSet = true;
  }
  if (command==fIexpoCmd) {
    fIexpo = fIexpoCmd->GetNewDoubleValue(newValue);
    fIexpoIsSet = true;
  }
  
  if (fZminIsSet && fZmaxIsSet && fIbasesIsSet && fIexpoIsSet) {
    if (!fICellCreated) {
      fITubeFactory->AddCell(fCellName, fZmin, fZmax);
      fICellCreated = true;
    }
    fITubeFactory->SetImportance(fCellName, fIbase, fIexpo);
  }

}










