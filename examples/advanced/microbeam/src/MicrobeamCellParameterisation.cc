// -------------------------------------------------------------------
// $Id: MicrobeamCellParameterisation.cc,v 1.1 2006-04-06 15:32:44 sincerti Exp $
// -------------------------------------------------------------------

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

#include "MicrobeamCellParameterisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MicrobeamCellParameterisation::MicrobeamCellParameterisation
(G4int NoBoxes, G4float DimBoxX, G4float DimBoxY, G4float DimBoxZ, 
 G4Material * nucleus,  G4Material * cytoplasm)
{
   nucleusMaterial = nucleus;
   cytoplasmMaterial = cytoplasm;
   
   NoCellBoxes = NoBoxes;
   DimCellBoxX = DimBoxX;
   DimCellBoxY = DimBoxY;
   DimCellBoxZ = DimBoxZ;
 
   mapCell  = new G4ThreeVector[NoCellBoxes];
   material = new G4float[NoCellBoxes];
 
   G4int ncols,nlines;
   G4int shiftX, shiftY, shiftZ;
   G4float x,y,z,mat,tmp,sizeZ;  
   ncols=0; nlines=0;
    
   // READ PHANTOM 
   FILE *fMap;
   fMap = fopen("phantom.dat","r");
   while (1) 
   {  
      if (nlines >=0 && nlines <=1) ncols = fscanf(fMap,"%f %f %f",&tmp,&tmp,&sizeZ);
      if (nlines ==2) ncols = fscanf(fMap,"%i %i %i",&shiftX,&shiftY,&shiftZ); // VOXEL SHIFT IN Z ASSUMED TO BE NEGATIVE
      if (nlines >2 ) ncols = fscanf(fMap,"%f %f %f %f %f",&x,&y,&z,&mat,&tmp);
      if (ncols < 0) break;

      G4ThreeVector v(x+shiftX,y+shiftY,z-1500/sizeZ-shiftZ); // VOXEL SHIFT TO CENTER PHANTOM
      
      if (nlines>2) 
      	{
	  mapCell[nlines-3]=v; 
	  material[nlines-3]=mat;
	}
      nlines++;    
   }
   fclose(fMap);
   
  // NUCLEUS IN GREEN 
  nucleusAttributes = new G4VisAttributes;
  nucleusAttributes->SetColour(G4Colour(0,1,0));
  nucleusAttributes->SetForceSolid(false);

  // CYTOPLASM IN RED
  cytoplasmAttributes = new G4VisAttributes;
  cytoplasmAttributes->SetColour(G4Colour(1,0,0));
  cytoplasmAttributes->SetForceSolid(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MicrobeamCellParameterisation::~MicrobeamCellParameterisation()
{
delete[] mapCell;
delete[] material;
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
    if( material[copyNo] == 2 ) 
        {
	 physVol->SetName("physicalNucleus");
	 physVol->GetLogicalVolume()->SetVisAttributes( nucleusAttributes );
	 return nucleusMaterial;
	} 
	
    else if( material[copyNo] == 1 )
	{
	 physVol->SetName("physicalCytoplasm");
	 physVol->GetLogicalVolume()->SetVisAttributes( cytoplasmAttributes );
	 return cytoplasmMaterial;
	} 

    return physVol->GetLogicalVolume()->GetMaterial();
}


