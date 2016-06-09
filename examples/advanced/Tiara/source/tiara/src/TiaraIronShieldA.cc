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
// $Id: TiaraIronShieldA.cc,v 1.3 2003/06/25 09:13:04 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
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
