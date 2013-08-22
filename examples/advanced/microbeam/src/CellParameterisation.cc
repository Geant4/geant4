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

#include "CellParameterisation.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CellParameterisation::CellParameterisation
(G4int NoBoxes, G4float DimBoxX, G4float DimBoxY, G4float DimBoxZ, 
 G4Material * nucleus1,  G4Material * cytoplasm1,
 G4Material * nucleus2,  G4Material * cytoplasm2,
 G4Material * nucleus3,  G4Material * cytoplasm3
 )
{
   fNucleusMaterial1   = nucleus1;
   fCytoplasmMaterial1 = cytoplasm1;
   fNucleusMaterial2   = nucleus2;
   fCytoplasmMaterial2 = cytoplasm2;
   fNucleusMaterial3   = nucleus3;
   fCytoplasmMaterial3 = cytoplasm3;

   fNoCellBoxes = NoBoxes;
   fDimCellBoxX = DimBoxX;
   fDimCellBoxY = DimBoxY;
   fDimCellBoxZ = DimBoxZ;
 
   fMapCell  = new G4ThreeVector[fNoCellBoxes];
   fMaterial = new G4float[fNoCellBoxes];
   fMass     = new G4float[fNoCellBoxes];
 
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

      G4ThreeVector v(x+shiftX,y+shiftY,z-1500/sizeZ-shiftZ); // VOXEL SHIFT IN ORDER TO CENTER PHANTOM
      
      if (nlines>4) 
      {
	  fMapCell[nlines-5]=v; 
	  fMaterial[nlines-5]=mat;
	  fMass[nlines-5]=den;
      }	  

      nlines++;    
   }
   fclose(fMap);
   
  // NUCLEUS IN GREEN 
  
  fNucleusAttributes1 = new G4VisAttributes;
  fNucleusAttributes1->SetColour(G4Colour(0,.8,0));
  fNucleusAttributes1->SetForceSolid(false);

  fNucleusAttributes2 = new G4VisAttributes;
  fNucleusAttributes2->SetColour(G4Colour(0,.9,0));
  fNucleusAttributes2->SetForceSolid(false);

  fNucleusAttributes3 = new G4VisAttributes;
  fNucleusAttributes3->SetColour(G4Colour(0,1,0));
  fNucleusAttributes3->SetForceSolid(false);

  // CYTOPLASM IN RED
  
  fCytoplasmAttributes1 = new G4VisAttributes;
  fCytoplasmAttributes1->SetColour(G4Colour(1,0,0));
  fCytoplasmAttributes1->SetForceSolid(false);

  fCytoplasmAttributes2 = new G4VisAttributes; // nucleoli in yellow
  fCytoplasmAttributes2->SetColour(G4Colour(1.,1.,0));
  fCytoplasmAttributes2->SetForceSolid(false);

  fCytoplasmAttributes3 = new G4VisAttributes;
  fCytoplasmAttributes3->SetColour(G4Colour(1,0,0));
  fCytoplasmAttributes3->SetForceSolid(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CellParameterisation::~CellParameterisation()
{
  delete[] fMapCell;
  delete[] fMaterial;
  delete[] fMass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CellParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  G4ThreeVector
    origin(
      fMapCell[copyNo].x()*fDimCellBoxX*2,
      fMapCell[copyNo].y()*fDimCellBoxY*2,
      fMapCell[copyNo].z()*fDimCellBoxZ*2);

  physVol->SetTranslation(origin);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CellParameterisation::ComputeDimensions
(G4Box& /*trackerChamber*/, const G4int /*copyNo*/, const G4VPhysicalVolume*) const
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material*
CellParameterisation::ComputeMaterial(const G4int copyNo,
                                               G4VPhysicalVolume* physVol,
                                               const G4VTouchable*)
{
    if( fMaterial[copyNo] == 2 ) // fMaterial 2 is nucleus
        {
	 if( fMass[copyNo] == 1 ) 
	 	{
	 	physVol->SetName("physicalNucleus");
	 	physVol->GetLogicalVolume()->SetVisAttributes( fNucleusAttributes1 );
		return fNucleusMaterial1;
		}
	 if( fMass[copyNo] == 2 ) 
	 	{
	 	physVol->SetName("physicalNucleus");
	 	physVol->GetLogicalVolume()->SetVisAttributes( fNucleusAttributes2 );
		return fNucleusMaterial2;
		}
	 if( fMass[copyNo] == 3 ) 
	 	{
	 	physVol->SetName("physicalNucleus");
	 	physVol->GetLogicalVolume()->SetVisAttributes( fNucleusAttributes3 );
		return fNucleusMaterial3;
		}
	} 
	
    else if( fMaterial[copyNo] == 1 ) // fMaterial 1 is cytoplasm
	{
	 if( fMass[copyNo] == 1 ) 
	 	{
	 	physVol->SetName("physicalCytoplasm");
	 	physVol->GetLogicalVolume()->SetVisAttributes( fCytoplasmAttributes1 );
	 	return fCytoplasmMaterial1;
		}
	 if( fMass[copyNo] == 2 ) 
	 	{
	 	physVol->SetName("physicalNucleus"); // nucleoli
	 	physVol->GetLogicalVolume()->SetVisAttributes( fCytoplasmAttributes2 );
	 	return fCytoplasmMaterial2;
		}
	 if( fMass[copyNo] == 3 ) 
	 	{
	 	physVol->SetName("physicalCytoplasm");
	 	physVol->GetLogicalVolume()->SetVisAttributes( fCytoplasmAttributes3 );
	 	return fCytoplasmMaterial3;
		}
	} 

    return physVol->GetLogicalVolume()->GetMaterial();
}
