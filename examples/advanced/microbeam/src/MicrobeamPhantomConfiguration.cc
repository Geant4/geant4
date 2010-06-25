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
// -------------------------------------------------------------------
// $Id: MicrobeamPhantomConfiguration.cc,v 1.9 2010-06-25 09:41:03 gunter Exp $
// -------------------------------------------------------------------

#include "MicrobeamPhantomConfiguration.hh"

G4int MicrobeamPhantomConfiguration::phantomTotalPixels = 0;
G4int MicrobeamPhantomConfiguration::nucleusTotalPixels = 0;
G4int MicrobeamPhantomConfiguration::cytoplasmTotalPixels = 0;
G4float MicrobeamPhantomConfiguration::dx = 0;
G4float MicrobeamPhantomConfiguration::dy = 0;
G4float MicrobeamPhantomConfiguration::dz = 0;
G4float MicrobeamPhantomConfiguration::nucleusMass = 0;
G4float MicrobeamPhantomConfiguration::cytoplasmMass = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MicrobeamPhantomConfiguration::MicrobeamPhantomConfiguration() {
Initialize();
}

MicrobeamPhantomConfiguration::~MicrobeamPhantomConfiguration()
{
  delete[] voxelThreeVector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int MicrobeamPhantomConfiguration::Initialize() {

  G4int ncols;
  G4float vx, vy, vz, tmp, density;
  G4int den, mat;
  G4float denCyto1, denCyto2, denCyto3, denNucl1, denNucl2, denNucl3;
  FILE* fMap;
  
  phantomTotalPixels=0;
  nucleusTotalPixels=0;
  cytoplasmTotalPixels=0;
  dx=0;
  dy=0;
  dz=0;
  nucleusMass=0;
  cytoplasmMass=0;
  density=0;

  // READ PHANTOM PARAMETERS
  fMap = fopen("phantom.dat","r");

  ncols = fscanf(fMap,"%i %i %i",&phantomTotalPixels, &nucleusTotalPixels, &cytoplasmTotalPixels);
  ncols = fscanf(fMap,"%f %f %f",&dx, &dy, &dz);
  ncols = fscanf(fMap,"%f %f %f",&tmp, &tmp, &tmp);
  ncols = fscanf(fMap,"%f %f %f",&denCyto1, &denCyto2, &denCyto3);
  ncols = fscanf(fMap,"%f %f %f",&denNucl1, &denNucl2, &denNucl3);
  dx = dx * micrometer;
  dy = dy * micrometer;
  dz = dz * micrometer;
  voxelThreeVector = new G4ThreeVector [phantomTotalPixels];

  for (G4int i=0; i<phantomTotalPixels; i++) 
  { 
    ncols = fscanf(fMap,"%f %f %f %i %i %f",&vx, &vy, &vz, &mat, &den, &tmp);

    if (std::abs(mat-2)<1.e-30) // NUCLEUS
    	{
	  if (std::abs(den-1)<1.e-30) density = denNucl1*(g/cm3);
	  if (std::abs(den-2)<1.e-30) density = denNucl2*(g/cm3);
	  if (std::abs(den-3)<1.e-30) density = denNucl3*(g/cm3);
	  nucleusMass   = nucleusMass   + density * dx * dy * dz ;
    	}

    if (std::abs(mat-1)<1.e-30) // CYTOPLASM
    	{ 
	  if (std::abs(den-1)<1e-30) density = denCyto1*(g/cm3);
	  if (std::abs(den-2)<1e-30) density = denCyto2*(g/cm3);
	  if (std::abs(den-3)<1e-30) density = denCyto3*(g/cm3);
	  cytoplasmMass = cytoplasmMass + density * dx * dy * dz ;
	}
    
    G4ThreeVector v(vx,vy,vz);
    voxelThreeVector[i] = v;
  }

  fclose(fMap);

  return 0;
}


