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
// $Id: RE02PhantomParameterisation.cc,v 1.1 2005/11/24 01:44:18 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
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






