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
// $Id: G4IWorldTubeMessenger.cc,v 1.2 2002-07-12 10:40:44 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4IWorldTubeMessenger.cc
//
// ----------------------------------------------------------------------

#include "G4IWorldTubeMessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4Tubs.hh"
#include "G4ImportanceGeometryConstructor.hh"
#include "G4LogicalVolume.hh"
#include "G4VIStore.hh"
#include "G4ITubeFactory.hh"

G4IWorldTubeMessenger::
G4IWorldTubeMessenger(G4ImportanceGeometryConstructor *igc):
  fIGconst(igc){
  fRadiusCmd = 
    new G4UIcmdWithADoubleAndUnit("/imp/worldvolume/tube/radius", this);
  fRadiusIsSet =false;
  fHalfHightCmd = 
    new G4UIcmdWithADoubleAndUnit("/imp/worldcolume/tube/halfwidth", this);
  fHalfhightIsSet = false;
  fITubeFactory = 0;
}

void G4IWorldTubeMessenger::
SetNewValue(G4UIcommand * command, G4String newValue){
  
  if (command==fRadiusCmd) {
    fRadius = fRadiusCmd->GetNewDoubleValue(newValue);
    fRadiusIsSet = true;
  }
  if (command==fHalfHightCmd) {
    fHalfHight = fHalfHightCmd->GetNewDoubleValue(newValue);
    fHalfhightIsSet = true;
  }
  
  if (fRadiusIsSet && fHalfhightIsSet) {
    ConstructWorldSolid();
    fIGconst->ConstructWorldVolume(fWorldSolid);
    ConstructICellFactory(fIGconst->GetLogicalWorld(),
			  fIGconst->GetIStore());
  }
  
}

void G4IWorldTubeMessenger::ConstructWorldSolid() {
  if (!fRadiusIsSet || !fHalfhightIsSet) {
    Error("!fRadiusIsSet || !fHalfhightIsSet");
  }
  if (!fWorldSolid) {
    G4String tubename("imp_WorldTube");
    fWorldSolid = new G4Tubs(tubename,
			     0.,
			     fRadius,
			     fHalfHight,
			     0.*deg,
			     360.*deg);
    G4cout << "G4IWorldTubeMessenger:: constructed World solid, with" 
	   << G4endl
	   <<"   radius    = " << fRadius << G4endl
	   <<"   halfhight = " 
	   << fHalfHight  << G4endl;
  }
}

void G4IWorldTubeMessenger::ConstructICellFactory(G4LogicalVolume *lv,
						  G4VIStore *is){
  if (!fRadiusIsSet || !fHalfhightIsSet) {
    Error("!fRadiusIsSet || !fHalfhightIsSet");
  }
  if (!fITubeFactory) {
    fITubeFactory = new G4ITubeFactory(lv,
				       is,
				       fRadius, 
				       fHalfHight);
  }
}

