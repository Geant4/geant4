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
// $Id: TiaraIronShieldA.cc,v 1.4 2006/06/29 15:45:06 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//

#include "TiaraIronShieldA.hh"
#include "TiaraDimensions.hh"
#include "TiaraMaterials.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "TiaraPart.hh"
#include "G4VisAttributes.hh"


TiaraIronShieldA::
TiaraIronShieldA(TiaraMaterials &mfac,
		 const TiaraDimensions &tiaraDimensions){
  G4double halfHight((tiaraDimensions.distTargetEndA - 
		      tiaraDimensions.distTargetWall) / 2);
  G4Tubs *ironAtub = new G4Tubs("ironAtub",
				   tiaraDimensions.pipeRadius,
				   tiaraDimensions.radiusIronA,
				   halfHight,
				   0,
				   360 * deg);
  G4LogicalVolume *logIronA = 
    new G4LogicalVolume(ironAtub, 
			mfac.GetMaterial("iron"),
			"IronA");

  G4VisAttributes* pIronAVis = new 
    G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  //  pIronAVis->SetForceSolid(true);
  logIronA->SetVisAttributes(pIronAVis);
  
  TiaraPart tiaraPart;
  tiaraPart.logVol = logIronA;
  tiaraPart.pos.setX(0);
  tiaraPart.pos.setY(0);
  tiaraPart.pos.setZ(tiaraDimensions.targetPosZ +  
		     tiaraDimensions.distTargetWall + halfHight);
  fTiaraParts.push_back(tiaraPart);
}

TiaraIronShieldA::~TiaraIronShieldA(){
}

TiaraParts TiaraIronShieldA::GetParts(){
  return fTiaraParts;
}
