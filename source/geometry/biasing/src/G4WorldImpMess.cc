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
// $Id: G4WorldImpMess.cc,v 1.2 2002-07-12 10:40:44 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4WorldImpMess.cc
//
// ----------------------------------------------------------------------

#include "G4WorldImpMess.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4ImportanceGeometryConstructor.hh"

G4WorldImpMess::G4WorldImpMess(G4ImportanceGeometryConstructor *igc) 
{
  fIGConst = igc;

  fIbaseCmd = new G4UIcmdWithADouble("/imp/worldvolume/ibase", this);
  fIbaseIsSet = false;
  fIexpoCmd = new G4UIcmdWithADouble("/imp/worldvolume/iexpo", this);
  fIexpoIsSet = false;

};


void G4WorldImpMess::
SetNewValue(G4UIcommand * command, G4String newValue){
  if (command==fIbaseCmd) {
    fIbase = fIbaseCmd->GetNewDoubleValue(newValue);
    fIbaseIsSet = true;
  }
  if (command==fIexpoCmd) {
    fIexpo = fIexpoCmd->GetNewDoubleValue(newValue);
    fIexpoIsSet = true;
  }
  if (fIbaseIsSet && fIexpoIsSet) {
    fIGConst->SetWorldImportance(fIbase, fIexpo);
  }

}
