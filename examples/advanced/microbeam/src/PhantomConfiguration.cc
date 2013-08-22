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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// If you use this example, please cite the following publication:
// Rad. Prot. Dos. 133 (2009) 2-11

#include "PhantomConfiguration.hh"
#include "G4SystemOfUnits.hh"

G4int PhantomConfiguration::fPhantomTotalPixels = 0;
G4int PhantomConfiguration::fNucleusTotalPixels = 0;
G4int PhantomConfiguration::fCytoplasmTotalPixels = 0;
G4float PhantomConfiguration::fDx = 0;
G4float PhantomConfiguration::fDy = 0;
G4float PhantomConfiguration::fDz = 0;
G4float PhantomConfiguration::fNucleusMass = 0;
G4float PhantomConfiguration::fCytoplasmMass = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhantomConfiguration::PhantomConfiguration() 
{
  Initialize();
}

PhantomConfiguration::~PhantomConfiguration()
{
  delete[] fVoxelThreeVector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int PhantomConfiguration::Initialize() 
{
  G4float vx, vy, vz, tmp, density;
  G4int den, mat;
  G4float denCyto1, denCyto2, denCyto3, denNucl1, denNucl2, denNucl3;
  
  FILE* fMap;
  
  density=0;

  fPhantomTotalPixels=0;
  fNucleusTotalPixels=0;
  fCytoplasmTotalPixels=0;
  fDx=0;
  fDy=0;
  fDz=0;
  fNucleusMass=0;
  fCytoplasmMass=0;
  
  // READ PHANTOM PARAMETERS
  fMap = fopen("phantom.dat","r");

  fscanf(fMap,"%i %i %i",&fPhantomTotalPixels, &fNucleusTotalPixels, &fCytoplasmTotalPixels);
  fscanf(fMap,"%f %f %f",&fDx, &fDy, &fDz);
  fscanf(fMap,"%f %f %f",&tmp, &tmp, &tmp);
  fscanf(fMap,"%f %f %f",&denCyto1, &denCyto2, &denCyto3);
  fscanf(fMap,"%f %f %f",&denNucl1, &denNucl2, &denNucl3);
  fDx = fDx * micrometer;
  fDy = fDy * micrometer;
  fDz = fDz * micrometer;
  
  fVoxelThreeVector = new G4ThreeVector [fPhantomTotalPixels];

  for (G4int i=0; i<fPhantomTotalPixels; i++) 
  { 
   fscanf(fMap,"%f %f %f %i %i %f",&vx, &vy, &vz, &mat, &den, &tmp);

    if (std::abs(mat-2)<1.e-30) // NUCLEUS
    	{
	  if (std::abs(den-1)<1.e-30) density = denNucl1*(g/cm3);
	  if (std::abs(den-2)<1.e-30) density = denNucl2*(g/cm3);
	  if (std::abs(den-3)<1.e-30) density = denNucl3*(g/cm3);
	  fNucleusMass   = fNucleusMass   + density * fDx * fDy * fDz ;
    	}

    if (std::abs(mat-1)<1.e-30) // CYTOPLASM
    	{ 
	  if (std::abs(den-1)<1e-30) density = denCyto1*(g/cm3);
	  if (std::abs(den-2)<1e-30) density = denCyto2*(g/cm3);
	  if (std::abs(den-3)<1e-30) density = denCyto3*(g/cm3);
	  fCytoplasmMass = fCytoplasmMass + density * fDx * fDy * fDz ;
	}
    
    G4ThreeVector v(vx,vy,vz);
    fVoxelThreeVector[i] = v;
  }

  fclose(fMap);

  return 0;
}

