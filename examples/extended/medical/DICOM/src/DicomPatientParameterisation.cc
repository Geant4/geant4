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

#include <string>

using namespace std;

DicomPatientParameterisation::~DicomPatientParameterisation()
{
  delete Attributes_Adipose;
  delete Attributes_LungEXhale;
  delete Attributes_Breast;
  delete Attributes_Phantom;
  delete Attributes_Muscle;
  delete Attributes_Liver;
  delete Attributes_TrabecularBone;
  delete Attributes_LungINhale;
  delete Attributes_DenseBone;
  delete Attributes_air;
}

DicomPatientParameterisation::DicomPatientParameterisation(G4int NoVoxels, double max_density, double min_density ,
							   G4Material* lunginhale,
							   G4Material* lungexhale,
							   G4Material* adipose_tissue,
							   G4Material* breast,
							   G4Material* phantom,
							   G4Material* muscle,
							   G4Material* liver,
							   G4Material* dense_bone,
							   G4Material* trabecular_bone)
{
  DicomConfiguration* ReadConfiguration = new DicomConfiguration;
  ReadConfiguration->ReadDataFile();					// images must have the same dimension
  G4int totalNumberOfFile = ReadConfiguration -> GetTotalNumberOfFile();
  for (int i=0;i<totalNumberOfFile;i++)
    {
      ReadConfiguration->ReadG4File( ReadConfiguration->GetListOfFile()[i] );
      MiddleLocationValue=MiddleLocationValue+ReadConfiguration->GetSliceLocation();
    }
  MiddleLocationValue=MiddleLocationValue/totalNumberOfFile;
 

    

  Attributes_LungINhale	= new G4VisAttributes;
  Attributes_LungINhale->SetColour(red=0.217/2,green=0.217/2,blue=0.217/2,alpha=1.);
  Attributes_LungINhale->SetForceSolid(true);

  Attributes_LungEXhale	= new G4VisAttributes;
  Attributes_LungEXhale->SetColour(red=0.508/2,green=0.508/2,blue=0.508/2,alpha=1.);
  Attributes_LungEXhale->SetForceSolid(true);

  Attributes_Adipose	= new G4VisAttributes;
  Attributes_Adipose->SetColour(red=0.967/2,green=0.967/2,blue=0.967/2,alpha=1.);
  Attributes_Adipose->SetForceSolid(true);

  Attributes_Breast	= new G4VisAttributes;
  Attributes_Breast->SetColour(red=0.99/2,green=0.99/2,blue=0.99/2,alpha=1.);
  Attributes_Breast->SetForceSolid(true);

  Attributes_Phantom	= new G4VisAttributes;
  Attributes_Phantom->SetColour(red=1.018/2,green=1.018/2,blue=1.018/2,alpha=1.);
  Attributes_Phantom->SetForceSolid(true);

  Attributes_Muscle	= new G4VisAttributes;
  Attributes_Muscle->SetColour(red=1.061/2,green=1.061/2,blue=1.061/2,alpha=1.);
  Attributes_Muscle->SetForceSolid(true);


  Attributes_Liver = new G4VisAttributes;
  Attributes_Liver->SetColour(red=1.071/2,green=1.071/2,blue=1.071/2,alpha=1.);
  Attributes_Liver->SetForceSolid(true);

  Attributes_TrabecularBone	= new G4VisAttributes;
  Attributes_TrabecularBone->SetColour(red=1.159/2,green=1.159/2,blue=1.159/2,alpha=1.);
  Attributes_TrabecularBone->SetForceSolid(true);

  Attributes_DenseBone	= new G4VisAttributes;
  Attributes_DenseBone->SetColour(red=1.575/2,green=1.575/2,blue=1.575/2,alpha=1.);
  Attributes_DenseBone->SetForceSolid(true);

  Attributes_air = new G4VisAttributes;
  Attributes_air->SetColour(red=0,green=0,blue=1,alpha=1.);
  Attributes_air->SetForceSolid(false);

  readData = fopen("Data.dat","r");
  fscanf(readData,"%s",compressionbuf);
  compression=atoi(compressionbuf);
  fscanf(readData,"%s",maxbuf);
  max=atoi(maxbuf);
  fscanf(readData,"%s",name);
  fclose(readData);

  sprintf(fullname,"%s.g4",name);
  readData =  fopen(fullname,"r");

  fscanf(readData,"%s %s",rowsbuf,columnsbuf);
  rows=atoi(rowsbuf);
  columns=atoi(columnsbuf);
  fscanf(readData,"%s %s",pixel_spacing_Xbuf,pixel_spacing_Ybuf);
  pixel_spacing_X=atof(pixel_spacing_Xbuf);
  pixel_spacing_Y=atof(pixel_spacing_Ybuf);
  fscanf(readData,"%s",SliceTicknessbuf);
  SliceTickness=atoi(SliceTicknessbuf);
  fscanf(readData,"%s",SliceLocationbuf);
  SliceLocation=atof(SliceLocationbuf);
  fscanf(readData,"%s",compressionbuf);
  compression=atoi(compressionbuf);
    
  GetDensity( max_density , min_density );

  P_lung_exhale=lungexhale;
  P_lung_inhale=lunginhale;
  P_adipose= adipose_tissue;
  P_breast=breast;
  P_phantom=phantom;
  P_muscle=muscle;
  P_liver=liver;
  P_dense_bone=dense_bone;
  P_trabecular_bone=trabecular_bone;
}

void DicomPatientParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{

  G4ThreeVector origin(PatientPlacementX[copyNo]*mm ,PatientPlacementY[copyNo]*mm ,PatientPlacementZ[copyNo]*mm-MiddleLocationValue*mm-SliceTickness/2*mm );
  physVol->SetTranslation(origin);
}


void DicomPatientParameterisation::ComputeDimensions(G4Box& Voxels, const G4int copyNo, const G4VPhysicalVolume* physVol) const
{
  Voxels.SetXHalfLength((pixel_spacing_X*compression/2.0)*mm);
  Voxels.SetYHalfLength((pixel_spacing_Y*compression/2.0)*mm);
  Voxels.SetZHalfLength((SliceTickness/2.0)*mm);
}

G4Material*  DicomPatientParameterisation::ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* physVol)
{
  bool detecting = true;
  if (detecting == false )
    {
      // OutOfBoundaries
    }
  else
    {
      if ( Density[copyNo] >= 0.207 && Density[copyNo] <= 0.227 )
        {
	  physVol->SetName("Physical_LungINhale");
	  physVol->GetLogicalVolume()->SetVisAttributes( Attributes_LungINhale );
	  return P_lung_inhale;
        }
      else if ( Density[copyNo] >= 0.481 && Density[copyNo] <= 0.534 )
        {
	  physVol->SetName("Physical_LungEXhale");
	  physVol->GetLogicalVolume()->SetVisAttributes( Attributes_LungEXhale );
	  return P_lung_exhale;
        }
      else if ( Density[copyNo] >= 0.919 && Density[copyNo] <= 0.979 )
        {
	  physVol->SetName("Physical_Adipose");
	  physVol->GetLogicalVolume()->SetVisAttributes( Attributes_Adipose );
	  return P_adipose;
        }
      else if ( Density[copyNo] > 0.979 && Density[copyNo] <= 1.004 )
        {
	  physVol->SetName("Physical_Breast");
	  physVol->GetLogicalVolume()->SetVisAttributes( Attributes_Breast );
	  return P_breast;
        }
      else if ( Density[copyNo] > 1.004 && Density[copyNo] <= 1.043 )
        {
	  physVol->SetName("Physical_Phantom");
	  physVol->GetLogicalVolume()->SetVisAttributes( Attributes_Phantom );
	  return P_phantom;
        }
      else if ( Density[copyNo] > 1.043 && Density[copyNo] <= 1.109 )
        {
	  physVol->SetName("Physical_Muscle");
	  physVol->GetLogicalVolume()->SetVisAttributes( Attributes_Muscle );
	  return P_muscle;
        }
      else if ( Density[copyNo] > 1.109 && Density[copyNo] <= 1.113 )
        {
	  physVol->SetName("Physical_Liver");
	  physVol->GetLogicalVolume()->SetVisAttributes( Attributes_Liver );
	  return P_liver;
        }
      else if ( Density[copyNo] > 1.113 && Density[copyNo] <= 1.217 )
        {
	  physVol->SetName("Physical_TrabecularBone");
	  physVol->GetLogicalVolume()->SetVisAttributes( Attributes_TrabecularBone );
	  return P_trabecular_bone;
        }
      else if ( Density[copyNo] > 1.496 && Density[copyNo] <= 1.654 )
        {
	  physVol->SetName("Physical_DenseBone");
	  physVol->GetLogicalVolume()->SetVisAttributes( Attributes_DenseBone );
	  return P_dense_bone;
        }

    }
  return physVol->GetLogicalVolume()->GetMaterial();
}

void DicomPatientParameterisation::GetDensity(double maxdensity , double mindensity)
{
  DicomConfiguration* ReadConfiguration = new DicomConfiguration;
  ReadConfiguration->ReadDataFile();

  int copy_counter = 0;
  G4int totalNumberOfFile = ReadConfiguration->GetTotalNumberOfFile();
  for (int z=0;z<totalNumberOfFile;z++)
    {


      ReadConfiguration->ReadG4File( ReadConfiguration->GetListOfFile()[z] );

      G4int compressionValue = ReadConfiguration->GetCompressionValue(); 
      lenr=abs(rows/compressionValue);
      lenc=abs(columns/compressionValue);

      int i=0;
      for (int j=1;j<=lenr;j++)//lenr
        {
	  for (int w=1;w<=lenc;w++)//lenc
            {
	      i++;
	      if ( ReadConfiguration->DensityValue[i] != -1 )
                {
		  if ( ReadConfiguration->DensityValue[i] >= mindensity && ReadConfiguration->DensityValue[i] <= maxdensity )
                    {
		      Density.push_back( ReadConfiguration->DensityValue[i] );
		      copy_counter++;
		      PatientPlacementX.push_back( (ReadConfiguration->IsCompressionUsed()*( ((lenc*ReadConfiguration->GetXPixelSpacing())/2)-(ReadConfiguration->GetYPixelSpacing()/2+(w-1)*ReadConfiguration->GetYPixelSpacing()) ) )*mm );
		      PatientPlacementY.push_back( (ReadConfiguration->IsCompressionUsed()*( ((lenr*ReadConfiguration->GetXPixelSpacing())/2)-(ReadConfiguration->GetYPixelSpacing()/2+(j-1)*ReadConfiguration->GetYPixelSpacing()) ) )*mm );
		      PatientPlacementZ.push_back( (ReadConfiguration->GetSliceLocation() + ReadConfiguration->GetSliceThickness()/2)*mm );
                    }
                }
            }
            
        }
    }

}

