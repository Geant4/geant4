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
// $Id$
// -------------------------------------------------------------------

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

#include "MicrobeamCellParameterisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MicrobeamCellParameterisation::MicrobeamCellParameterisation
(G4int NoBoxes, G4float DimBoxX, G4float DimBoxY, G4float DimBoxZ, 
 G4Material * nucleus1,  G4Material * cytoplasm1,
 G4Material * nucleus2,  G4Material * cytoplasm2,
 G4Material * nucleus3,  G4Material * cytoplasm3
 )
{
   nucleusMaterial1 = nucleus1;
   cytoplasmMaterial1 = cytoplasm1;
   nucleusMaterial2 = nucleus2;
   cytoplasmMaterial2 = cytoplasm2;
   nucleusMaterial3 = nucleus3;
   cytoplasmMaterial3 = cytoplasm3;

   NoCellBoxes = NoBoxes;
   DimCellBoxX = DimBoxX;
   DimCellBoxY = DimBoxY;
   DimCellBoxZ = DimBoxZ;
 
   mapCell  = new G4ThreeVector[NoCellBoxes];
   material = new G4float[NoCellBoxes];
   mass     = new G4float[NoCellBoxes];
 
   G4int ncols,nlines;
   G4int shiftX, shiftY, shiftZ;
   G4float x,y,z,mat,den,tmp,sizeZ;  
   ncols=0; nlines=0;
    
   // READ PHANTOM 
   FILE *fMap;
   fMap = fopen("phantom.dat","r");
   while (1) 
   {  
      if (nlines >= 0  && nlines <=1 ) ncols = fscanf(fMap,"%f %f %f",&tmp,&tmp,&sizeZ);
      if (nlines == 2) ncols = fscanf(fMap,"%i %i %i",&shiftX,&shiftY,&shiftZ); // VOXEL SHIFT IN Z ASSUMED TO BE NEGATIVE
      if (nlines == 3) ncols = fscanf(fMap,"%f %f %f",&tmp,&tmp,&tmp);
      if (nlines == 4) ncols = fscanf(fMap,"%f %f %f",&tmp,&tmp,&tmp);
      if (nlines >  4) ncols = fscanf(fMap,"%f %f %f %f %f %f",&x,&y,&z,&mat,&den,&tmp);
      if (ncols  <  0) break;

      G4ThreeVector v(x+shiftX,y+shiftY,z-1500/sizeZ-shiftZ); // VOXEL SHIFT TO CENTER PHANTOM
      
      if (nlines>4) 
      {
	  mapCell[nlines-5]=v; 
	  material[nlines-5]=mat;
	  mass[nlines-5]=den;
      }	  

      nlines++;    
   }
   fclose(fMap);
   
  // NUCLEUS IN GREEN 
  nucleusAttributes1 = new G4VisAttributes;
  nucleusAttributes1->SetColour(G4Colour(0,.8,0));
  nucleusAttributes1->SetForceSolid(false);

  nucleusAttributes2 = new G4VisAttributes;
  nucleusAttributes2->SetColour(G4Colour(0,.9,0));
  nucleusAttributes2->SetForceSolid(false);

  nucleusAttributes3 = new G4VisAttributes;
  nucleusAttributes3->SetColour(G4Colour(0,1,0));
  nucleusAttributes3->SetForceSolid(false);

// CYTOPLASM IN RED
  cytoplasmAttributes1 = new G4VisAttributes;
  cytoplasmAttributes1->SetColour(G4Colour(1,0,0));
  cytoplasmAttributes1->SetForceSolid(false);

  cytoplasmAttributes2 = new G4VisAttributes; // nucleoli in yellow
  cytoplasmAttributes2->SetColour(G4Colour(1.,1.,0));
  cytoplasmAttributes2->SetForceSolid(false);

  cytoplasmAttributes3 = new G4VisAttributes;
  cytoplasmAttributes3->SetColour(G4Colour(1,0,0));
  cytoplasmAttributes3->SetForceSolid(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MicrobeamCellParameterisation::~MicrobeamCellParameterisation()
{
delete[] mapCell;
delete[] material;
delete[] mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MicrobeamCellParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  G4ThreeVector
    origin(
      mapCell[copyNo].x()*DimCellBoxX*2,
      mapCell[copyNo].y()*DimCellBoxY*2,
      mapCell[copyNo].z()*DimCellBoxZ*2);

  physVol->SetTranslation(origin);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MicrobeamCellParameterisation::ComputeDimensions
(G4Box& /*trackerChamber*/, const G4int /*copyNo*/, const G4VPhysicalVolume*) const
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material*
MicrobeamCellParameterisation::ComputeMaterial(const G4int copyNo,
                                               G4VPhysicalVolume* physVol,
                                               const G4VTouchable*)
{
    if( material[copyNo] == 2 ) // material 2 is nucleus
        {
	 if( mass[copyNo] == 1 ) 
	 	{
	 	physVol->SetName("physicalNucleus");
	 	physVol->GetLogicalVolume()->SetVisAttributes( nucleusAttributes1 );
		return nucleusMaterial1;
		}
	 if( mass[copyNo] == 2 ) 
	 	{
	 	physVol->SetName("physicalNucleus");
	 	physVol->GetLogicalVolume()->SetVisAttributes( nucleusAttributes2 );
		return nucleusMaterial2;
		}
	 if( mass[copyNo] == 3 ) 
	 	{
	 	physVol->SetName("physicalNucleus");
	 	physVol->GetLogicalVolume()->SetVisAttributes( nucleusAttributes3 );
		return nucleusMaterial3;
		}
	} 
	
    else if( material[copyNo] == 1 ) // material 1 is cytoplasm
	{
	 if( mass[copyNo] == 1 ) 
	 	{
	 	physVol->SetName("physicalCytoplasm");
	 	physVol->GetLogicalVolume()->SetVisAttributes( cytoplasmAttributes1 );
	 	return cytoplasmMaterial1;
		}
	 if( mass[copyNo] == 2 ) 
	 	{
	 	physVol->SetName("physicalNucleus"); // nucleoli
	 	physVol->GetLogicalVolume()->SetVisAttributes( cytoplasmAttributes2 );
	 	return cytoplasmMaterial2;
		}
	 if( mass[copyNo] == 3 ) 
	 	{
	 	physVol->SetName("physicalCytoplasm");
	 	physVol->GetLogicalVolume()->SetVisAttributes( cytoplasmAttributes3 );
	 	return cytoplasmMaterial3;
		}
	} 

    return physVol->GetLogicalVolume()->GetMaterial();
}


