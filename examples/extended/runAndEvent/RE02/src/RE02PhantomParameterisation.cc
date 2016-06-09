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
// $Id: RE02PhantomParameterisation.cc,v 1.2 2006/06/29 17:45:26 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//

#include "RE02PhantomParameterisation.hh"
#include "G4Box.hh"
#include "G4VPhysicalVolume.hh"

RE02PhantomParameterisation::RE02PhantomParameterisation(
		   const G4ThreeVector& motherSize,
		   const G4int nx, const G4int ny, const G4int nz)
    :G4VPVParameterisation(),fDxyzMother(motherSize),
     fNx(nx),fNy(ny),fNz(nz)
{
  // Voxel Size
  fDxyz.setX(fDxyzMother.x()/(G4double)fNx);
  fDxyz.setY(fDxyzMother.y()/(G4double)fNy);
  fDxyz.setZ(fDxyzMother.z()/(G4double)fNz);
  // Calculation of each segmented position, and fill it in vector.
  G4double offsetX = -fDxyzMother.x()+fDxyz.x();
  G4double offsetY = -fDxyzMother.y()+fDxyz.y();
  G4double offsetZ = -fDxyzMother.z()+fDxyz.z();
  for ( G4int ix = 0; ix < fNx; ix++){
      for ( G4int iy = 0; iy < fNy; iy++){
	  for ( G4int iz = 0; iz < fNz; iz++){
	      G4double x = offsetX+((G4double)ix)*fDxyz.x()*2.;
	      G4double y = offsetY+((G4double)iy)*fDxyz.y()*2.;
	      G4double z = offsetZ+((G4double)iz)*fDxyz.z()*2.;
	      G4ThreeVector pos(x,y,z);
	      fPositions.push_back(pos);
	  }
      }
  }

}

RE02PhantomParameterisation::~RE02PhantomParameterisation() {
}

void RE02PhantomParameterisation::ComputeTransformation(const G4int copyNo, 
				  G4VPhysicalVolume* physVol) const
{
    physVol->SetTranslation(fPositions[copyNo]);
}






