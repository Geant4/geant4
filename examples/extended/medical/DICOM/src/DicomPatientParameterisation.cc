//   $tigre.2@sympatico.ca, louis.archambault@phy.ulaval.ca
//   01/10/02

//*******************************************************
//
// DicomPatientParameterisation.cc :
//	- Parameterisation class for DICM images.
// 	- Placement and resizing of G4Box to do the
//	  final 3d geometry (DICM images)
//	- Setting good materials for each one of the
//	  G4Box
//	- In the future we hope to have multiple size
//	  G4Box to make the program faster, but for now
//	  all of them are of the same size ( there's alot
//	  of them !)
// 	- Definitions are in DicomPatientParameterisation.hh
//
// The code was written by :
//	Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
// For more information contact :
//	Louis Archambault louis.archambault@phy.ulaval.ca
// or
//	Luc Beaulieu beaulieu@phy.ulaval.ca
//
// Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//*******************************************************

#include "DicomPrimaryGeneratorAction.hh"
#include "DicomPatientParameterisation.hh"
#include "DicomConfiguration.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

DicomPatientParameterisation::DicomPatientParameterisation(G4int, // NoVoxels, 
							   G4double maxDensity, 
							   G4double minDensity ,
							   G4Material* lunginhale,
							   G4Material* lungexhale,
							   G4Material* adipose,
							   G4Material* breast,
							   G4Material* phantom,
							   G4Material* muscle,
							   G4Material* liver,
							   G4Material* denseBone,
							   G4Material* trabecularBone)
{ 
  lungExhale = lungexhale;
  lungInhale = lunginhale;
  adiposeTissue = adipose;
  breastTissue = breast;
  phantomTissue = phantom;
  muscleTissue = muscle;
  liverTissue = liver;
  denseBoneTissue = denseBone;
  trabecularBoneTissue = trabecularBone;
 
  DicomConfiguration* dicomConfiguration = new DicomConfiguration;
  dicomConfiguration->ReadDataFile();					
  // images must have the same dimension
  G4int totalNumberOfFile = dicomConfiguration -> GetTotalNumberOfFile();
 
  for ( G4int i = 0; i < totalNumberOfFile; i++)
    {
      dicomConfiguration->ReadG4File( dicomConfiguration->GetListOfFile()[i] );
      G4double sliceLocation = dicomConfiguration->GetSliceLocation();
      middleLocationValue = middleLocationValue + sliceLocation;
    }
  
  delete dicomConfiguration;

  middleLocationValue = middleLocationValue/totalNumberOfFile;   
  
  G4double red;
  G4double green;
  G4double blue;
  G4double alpha;
  
  attributeLungINhale = new G4VisAttributes;
  attributeLungINhale->SetColour(red=0.217/2,green=0.217/2,blue=0.217/2,alpha=1.);
  attributeLungINhale->SetForceSolid(true);

  attributeLungEXhale = new G4VisAttributes;
  attributeLungEXhale->SetColour(red=0.508/2,green=0.508/2,blue=0.508/2,alpha=1.);
  attributeLungEXhale->SetForceSolid(true);

  attributeAdipose = new G4VisAttributes;
  attributeAdipose->SetColour(red=0.967/2,green=0.967/2,blue=0.967/2,alpha=1.);
  attributeAdipose->SetForceSolid(true);

  attributeBreast = new G4VisAttributes;
  attributeBreast->SetColour(red=0.99/2,green=0.99/2,blue=0.99/2,alpha=1.);
  attributeBreast->SetForceSolid(true);

  attributePhantom = new G4VisAttributes;
  attributePhantom->SetColour(red=1.018/2,green=1.018/2,blue=1.018/2,alpha=1.);
  attributePhantom->SetForceSolid(true);

  attributeMuscle = new G4VisAttributes;
  attributeMuscle->SetColour(red=1.061/2,green=1.061/2,blue=1.061/2,alpha=1.);
  attributeMuscle->SetForceSolid(true);


  attributeLiver = new G4VisAttributes;
  attributeLiver->SetColour(red=1.071/2,green=1.071/2,blue=1.071/2,alpha=1.);
  attributeLiver->SetForceSolid(true);

  attributeTrabecularBone = new G4VisAttributes;
  attributeTrabecularBone->SetColour(red=1.159/2,green=1.159/2,blue=1.159/2,alpha=1.);
  attributeTrabecularBone->SetForceSolid(true);

  attributeDenseBone = new G4VisAttributes;
  attributeDenseBone->SetColour(red=1.575/2,green=1.575/2,blue=1.575/2,alpha=1.);
  attributeDenseBone->SetForceSolid(true);

  attributeAir = new G4VisAttributes;
  attributeAir->SetColour(red=0,green=0,blue=1,alpha=1.);
  attributeAir->SetForceSolid(false);
  
  char name[300];
  char maxBuf[300];
  char rowsBuf[300];
  char columnsBuf[300];
  char pixelSpacingXBuf[300];
  char pixelSpacingYBuf[300];
  char sliceThicknessBuf[300];  
  char sliceLocationBuf[300];
  char compressionBuf[300];
  char fullFileName[300];
  readData = G4std::fopen("Data.dat","r");
  G4std::fscanf(readData,"%s",compressionBuf);
  compression = atoi(compressionBuf);
  G4std::fscanf(readData,"%s",maxBuf);
  max = atoi(maxBuf);
  G4std::fscanf(readData,"%s",name);
  G4std::fclose(readData);

  G4std::sprintf(fullFileName,"%s.g4",name);
  readData = G4std::fopen(fullFileName,"r");

  G4std::fscanf(readData,"%s %s",rowsBuf,columnsBuf);
  rows = atoi(rowsBuf);
  columns = atoi(columnsBuf);
  G4std::fscanf(readData,"%s %s",pixelSpacingXBuf,pixelSpacingYBuf);
  pixelSpacingX = atof(pixelSpacingXBuf);
  pixelSpacingY = atof(pixelSpacingYBuf);
  G4std::fscanf(readData,"%s",sliceThicknessBuf);
  sliceThickness = atoi(sliceThicknessBuf);
  G4std::fscanf(readData,"%s",sliceLocationBuf);
  sliceLocation = atof(sliceLocationBuf);
  G4std::fscanf(readData,"%s",compressionBuf);
  compression = atoi(compressionBuf);
    
  GetDensity( maxDensity , minDensity );

 
}

DicomPatientParameterisation::~DicomPatientParameterisation()
{
  // visualisation attributes ...
  delete attributeAdipose;
  delete attributeLungEXhale;
  delete attributeBreast;
  delete attributePhantom;
  delete attributeMuscle;
  delete attributeLiver;
  delete attributeTrabecularBone;
  delete attributeLungINhale;
  delete attributeDenseBone;
  delete attributeAir;
  
  // materials ...
  delete trabecularBoneTissue;
  delete denseBoneTissue;
  delete liverTissue;
  delete muscleTissue;
  delete phantomTissue;
  delete breastTissue; 
  delete adiposeTissue;
  delete lungInhale;
  delete lungExhale;
}

void DicomPatientParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  G4double originZ = patientPlacementZ[copyNo]*mm-middleLocationValue*mm-sliceThickness/2;
  G4ThreeVector origin( patientPlacementX[copyNo]*mm, 
                        patientPlacementY[copyNo]*mm, 
                        originZ*mm );

  physVol->SetTranslation(origin);
}

void DicomPatientParameterisation::ComputeDimensions(G4Box& voxels, const G4int, const G4VPhysicalVolume*) const
{
  voxels.SetXHalfLength((pixelSpacingX * compression/2.0) * mm);
  voxels.SetYHalfLength((pixelSpacingY * compression/2.0) * mm);
  voxels.SetZHalfLength((sliceThickness / 2.0) * mm);
}

G4Material*  DicomPatientParameterisation::ComputeMaterial(const G4int copyNo,G4VPhysicalVolume* physVol)
{
  if ( density[copyNo] >= 0.207 && density[copyNo] <= 0.227 )
    {
      physVol->SetName("Physical_LungINhale");
      physVol->GetLogicalVolume()->SetVisAttributes( attributeLungINhale );
      return lungInhale;
    }
  else if ( density[copyNo] >= 0.481 && density[copyNo] <= 0.534 )
    {
      physVol->SetName("Physical_LungEXhale");
      physVol->GetLogicalVolume()->SetVisAttributes( attributeLungEXhale );
      return lungExhale;
    }
  else if ( density[copyNo] >= 0.919 && density[copyNo] <= 0.979 )
    {
      physVol->SetName("Physical_Adipose");
      physVol->GetLogicalVolume()->SetVisAttributes( attributeAdipose );
      return adiposeTissue;
    }
  else if ( density[copyNo] > 0.979 && density[copyNo] <= 1.004 )
    {
      physVol->SetName("Physical_Breast");
      physVol->GetLogicalVolume()->SetVisAttributes( attributeBreast );
      return breastTissue;
    }
  else if ( density[copyNo] > 1.004 && density[copyNo] <= 1.043 )
    {
      physVol->SetName("Physical_Phantom");
      physVol->GetLogicalVolume()->SetVisAttributes( attributePhantom );
      return phantomTissue;
    }
  else if ( density[copyNo] > 1.043 && density[copyNo] <= 1.109 )
    {
      physVol->SetName("Physical_Muscle");
      physVol->GetLogicalVolume()->SetVisAttributes( attributeMuscle );
      return muscleTissue;
    }
  else if ( density[copyNo] > 1.109 && density[copyNo] <= 1.113 )
    {
      physVol->SetName("Physical_Liver");
      physVol->GetLogicalVolume()->SetVisAttributes( attributeLiver );
      return liverTissue;
    }
  else if ( density[copyNo] > 1.113 && density[copyNo] <= 1.217 )
    {
      physVol->SetName("Physical_TrabecularBone");
      physVol->GetLogicalVolume()->SetVisAttributes( attributeTrabecularBone );
      return trabecularBoneTissue;
    }
  else if ( density[copyNo] > 1.496 && density[copyNo] <= 1.654 )
    {
      physVol->SetName("Physical_DenseBone");
      physVol->GetLogicalVolume()->SetVisAttributes( attributeDenseBone );
      return denseBoneTissue;
    }

  return physVol->GetLogicalVolume()->GetMaterial();
}

void DicomPatientParameterisation::GetDensity(G4double maxdensity, G4double mindensity)
{
  DicomConfiguration* dicomConfiguration = new DicomConfiguration;
  dicomConfiguration->ReadDataFile();

  G4int copyCounter = 0;
  G4int totalNumberOfFile = dicomConfiguration->GetTotalNumberOfFile();
  for ( G4int z = 0; z < totalNumberOfFile; z++ )
    {
      dicomConfiguration->ReadG4File( dicomConfiguration->GetListOfFile()[z] );
      G4int compressionValue = dicomConfiguration->GetCompressionValue(); 
      G4int lenRows = abs(rows/compressionValue);
      G4int lenColumns=abs(columns/compressionValue);

      G4int i = 0;
      for ( G4int j = 1; j <= lenRows; j++ )
        {
	  for ( G4int w = 1; w <= lenColumns; w++ )
            {
	      i++;
              G4double tissueDensity = dicomConfiguration->GetDensityValue(i);
	      if ( tissueDensity != -1 )
                {
		  if ( tissueDensity >= mindensity && tissueDensity <= maxdensity )
                    {
		      density.push_back( tissueDensity );
		      copyCounter++;
                      G4int isCompressionUsed = dicomConfiguration->IsCompressionUsed();
                      G4double xPixelSpacing =  dicomConfiguration->GetXPixelSpacing();
                      G4double yPixelSpacing =  dicomConfiguration->GetYPixelSpacing();
		      G4double slicePosition = dicomConfiguration->GetSliceLocation();      
                      G4double sliceThick =  dicomConfiguration->GetSliceThickness();
                      G4double xDimension = (lenColumns*xPixelSpacing)/2;
                      G4double yPixel = (yPixelSpacing/2+(w-1)*yPixelSpacing);
                      G4double yDimension = ((lenRows*xPixelSpacing)/2)-(yPixelSpacing/2+(j-1)*yPixelSpacing);
                      
                      patientPlacementX.push_back( ( isCompressionUsed*(xDimension- yPixel ) ) *mm );
                      patientPlacementY.push_back( ( isCompressionUsed* yDimension  ) *mm );
		      patientPlacementZ.push_back( ( slicePosition + sliceThick/2 ) *mm );
                    }
                }
            }            
        }
    }
  delete dicomConfiguration;
}

