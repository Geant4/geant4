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
// $Id: G4ITubeFactory.cc,v 1.2 2002-07-12 10:40:44 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ITubeFactory.cc
//
// ----------------------------------------------------------------------

#include "G4ITubeFactory.hh"
#include "G4ITubeMessenger.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VIStore.hh"
#include "G4Material.hh"

G4ITubeFactory::G4ITubeFactory(G4LogicalVolume *wl, G4VIStore *is,
				 G4double Radius, G4double HalfHight) :
  fWorldLogic(wl),
  fIStore(is),
  fRadius(Radius),
  fHalfHight(HalfHight),
  fITubeMessenger(new G4ITubeMessenger(this)),
  fGalactic(new G4Material("Galactic", 
			   1., 
			   1.01*g/mole, 
			   universe_mean_density,
			   kStateGas,
			   2.73*kelvin,
			   3.e-18*pascal))
{} 

void G4ITubeFactory::AddCell(const G4String &cellname, 
			      G4double zmin, G4double zmax){

  G4MapNameTube::iterator it = fMapNameTube.find(cellname);
  if (it != fMapNameTube.end()) {
    Error("cell: " + cellname + ", already exiosts");
  }
  
  G4double zwidthhalf = (zmax-zmin)/2;

  if (zwidthhalf<=0.) {
    Error("(zmax-zmin)/2 < 0" );
  }

  fMapNameTube[cellname] = new G4Tubs(cellname,
				     0.,
				     fRadius,
				     zwidthhalf,
				     0.*deg,
				     360.*deg);
  
  G4cout << "G4ITubeFactory:: constructed tubecell with: " << G4endl 
	 << "   zmin = " << zmin << G4endl 
	 << "   zmax = " <<  zmax << G4endl;


  fMapNameLogic[cellname] = 
    new G4LogicalVolume(fMapNameTube[cellname], fGalactic, 
                        cellname);
  
  fMapNamePhysical[cellname] = 
    new G4PVPlacement(0, G4ThreeVector(0,0,zmin+zwidthhalf), 
		      fMapNameLogic[cellname],
		      cellname, fWorldLogic, false, 0);
  
  fIStore->AddImportanceRegion(1., *(fMapNamePhysical[cellname]));
}

void G4ITubeFactory::SetImportance(const G4String &cellname, 
				    G4double b, G4double e){
  G4MapNamePhysical::iterator it = fMapNamePhysical.find(cellname);
  if (it==fMapNamePhysical.end() ){
    Error("cell: " + cellname + ", not defined");
  }
  fIStore->ChangeImportance(pow(b,e), *(it->second));
}

