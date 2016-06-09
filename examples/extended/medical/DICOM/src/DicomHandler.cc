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
// The code was written by :
//	*Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Université Laval, Québec (QC) Canada
//*******************************************************
//
//*******************************************************
//
//*******************************************************
//
// DicomHandler.cc :
//	- Handling of DICM images
// 	- Reading headers and pixels
//	- Transforming pixel to density and creating *.g4
//	  files
// 	- Definitions are in DicomHandler.hh
//*******************************************************

#include "DicomHandler.hh"
#include "globals.hh"
#include "G4ios.hh"
#include <strstream>
#include <streambuf>
#include "G4strstreambuf.hh"
#include <fstream>

DicomHandler::DicomHandler()
{
  compression = 0;
  max = 0;
}

G4int DicomHandler::readHeader(FILE *dicom, char filename2[300])
{
  G4int returnvalue = 0;
  char buffer[196];
  char  pixelSpacing[300];
  std::fread( buffer, 1, 128, dicom ); // The first 128 bytes 
                                         //are not important
  // Reads the "DICOM" letters
  std::fread( buffer, 1, 4, dicom );
  G4int readGroupId; //identify the kind of input data 
  G4int readElementId;//identify a particular type information

  // the elementLength say if a particular information (associated with a
  //given readGroupId+readElementId) is to be read in 2 or 4 bits
  G4int elementLength;
  G4int elementLength2;
  G4int elementLength3;
  char value[10000][300];
  // Read information up to the pixel data
  // note: it should be a while instead of a for
  for ( G4int i = 0; i <= 100000000; i++ )
    {
      //Reading groups and elements :
      std::fread(&readGroupId,1,2,dicom);
      std::fread(&readElementId,1,2,dicom);

      if (readGroupId == 0x7FE0) // beginning of the pixels
        {
	  std::fread( buffer, 1, 2,dicom); // Skip 2 reserved bytes
	  break;
        }

      std::fread(&elementLength,1,2,dicom);
      G4int tagDictionnary;
      G4int bitStored = 0;
      // If value representation (VR) is OB, OW, SQ, UN, 
      //the next length is 32 bits
      if ( elementLength == 16975 || 
           elementLength == 22351 || 
           elementLength == 20819 || 
           elementLength == 20053)
        {
	  //skip 2 reserved "bytes"
	  std::fread( buffer, 1, 2,dicom); // Skip 2 reserved bytes
	  std::fread(&elementLength3, 4, 1, dicom);
          // Reading length of the information
	  std::fread(&value[i],elementLength3,1,dicom);
          // Reading the information with
	  // (BIG) buffer : "value"
	  // Creating a tag to be identified afterward
	  tagDictionnary = readGroupId*0x10000 + readElementId;
        }
      else  // lenght is 16 bits :
        {
	  std::fread(&elementLength2,1,2,dicom);
	  std::fread(&value[i],elementLength2,1,dicom);
	  tagDictionnary = readGroupId*0x10000 + readElementId;
        }

      if (tagDictionnary == 0x00280010 ) // Number of Rows
        {
	  rows = *(G4int*)&value[i];
	  std::printf("[0x00280010] Rows -> %i\n",rows);
        }
      if (tagDictionnary == 0x00280011 ) // Number of columns
        {
	  columns = *(G4int*)&value[i];
	  std::printf("[0x00280011] Columns -> %i\n",columns);
        }
      if (tagDictionnary == 0x00280102 ) // High bits  ( not used )
        {
	  G4int highBits = *(G4int*)&value[i];
	  std::printf("[0x00280102] High bits -> %i\n",highBits);
        }
      if (tagDictionnary == 0x00280100 )  // Bits allocated ( not used )
        {
	  bitAllocated = *(G4int*)&value[i];
	  std::printf("[0x00280100] Bits allocated -> %i\n",bitAllocated);
	  bitAllocated = (bitAllocated)/8;
        }
      if (tagDictionnary == 0x00280101 )  //  Bits stored ( not used )
        {
	  bitStored = *(G4int*)&value[i];
	  std::printf("[0x00280101] Bits stord -> %i\n",bitStored);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00281053)  //  Rescale slope ( not used )
        { 
	  G4int rescaleSlope =  atoi( value[i] );
	  std::printf("[0x00281053] Rescale Slope -> %i\n",rescaleSlope);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00281052 )  // Rescalse intercept ( not used )
        {
	  G4int rescaleIntercept = atoi( value[i] );
	  std::printf("[0x00281052] Rescale Intercept -> %i\n", rescaleIntercept );
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00280103 )
        {
	  //  Pixel representation ( functions not design to read signed bits )
	  std::printf("[0x00280103] Pixel Representation -> %i\n",  atoi( value[i] ) );
	  if ( atoi(value[i]) == 1 )
            {
	      std::printf("### PIXEL REPRESENTATION = 1, BITS ARE SIGNED, ");
	      std::printf("DICOM READING SCAN FOR UNSIGNED VALUE, POSSIBLE ");
	      std::printf("ERROR !!!!!! -> \n");
            }
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00080008 ) //  Image type ( not used )
        {
	  std::printf("[0x00080008] Image Types -> %s\n",value[i]);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00283000 )   //  Modality LUT Sequence ( not used )
        {
	  std::printf("[0x00283000] Modality LUT Sequence SQ 1 -> %s\n",value[i]);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00283002 )  // LUT Descriptor ( not used )
        {
	  std::printf("[0x00283002] LUT Descriptor US or SS 3 -> %s\n",value[i]);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00283003 )  // LUT Explanation ( not used )
        {
	  std::printf("[0x00283003] LUT Explanation LO 1 -> %s\n",value[i]);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00283004 )  // Modality LUT ( not used )
        {
	  std::printf("[0x00283004] Modality LUT Type LO 1 -> %s\n",value[i]);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00283006 )  // LUT Data ( not used )
        {
	  std::printf("[0x00283006] LUT Data US or SS -> %s\n",value[i]);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00283010 )  // VOI LUT ( not used )
        {
	  std::printf("[0x00283010] VOI LUT Sequence SQ 1 -> %s\n",value[i]);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00280120 )  // Pixel Padding Value ( not used )
        {
	  std::printf("[0x00280120] Pixel Padding Value US or SS 1 -> %s\n",value[i]);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00280030 )  // Pixel Spacing
        {
	  std::printf("[0x00280030] Pixel Spacing (mm) -> %s\n",value[i]);
	  std::printf(pixelSpacing,"%s",value[i]);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00200037 )  // Image Orientation ( not used )
        {
	  std::printf("[0x00200037] Image Orientation (Patient) -> %s\n",value[i]);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00200032 )  // Image Position ( not used )
        {
	  std::printf("[0x00200032] Image Position (Patient,mm) -> %s\n",value[i]);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00180050 )  // Slice Tickness
        {
	  std::printf("[0x00180050] Slice Tickness (mm) -> %s\n",value[i]);
	  std::sprintf(sliceThickness,"%s",value[i]);//sliceThickness=value[i];
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00201041 )  // Slice Location
        {
	  std::printf("[0x00201041] Slice Location -> %s\n",value[i]);
	  sliceLocation = atof(value[i]);
	  bitStored = (bitStored)/8;
        }
      if (tagDictionnary == 0x00280004 ) 
      // Photometric Interpretation ( not used )
        {
	  std::printf("[0x00280004] Photometric Interpretation -> %s\n",value[i]);
	  bitStored = (bitStored)/8;
        }
    }

  // Creating files to store information
  char compressionbuf[100],maxbuf[100];
  char filename[300];
  compression = 0;
  max = 0;
  FILE* configuration;

  configuration = std::fopen("Data.dat","r");
  if ( configuration != 0 )
    {
      std::fscanf(configuration,"%s",compressionbuf);
      compression = atoi(compressionbuf);
      std::fscanf(configuration,"%s",maxbuf);
      max = atoi(maxbuf);
      std::fclose(configuration);
    }
  else
    {
      std::printf("### WARNING, file Data.dat not here !!!\n");
      exit(1);
    }
  FILE* data;
  std::sprintf(filename,"%s.dat",filename2);
  data = std::fopen(filename,"w+");
  // Note: the .dat files contain basic information on the images.

  char exception = '\\';
  G4bool toggle = false;
  G4int z = 0;
 
  for ( G4int y = 0; y <= 300; y++ )
    {
      if ( pixelSpacing[y] != exception )
        {
	  if (toggle == false)
	    pixelSpacingX[y] = pixelSpacing[y];
	  if (toggle == true)
            {
	      pixelSpacingY[z] = pixelSpacing[y];
	      z++;
            }
        }
      else if ( pixelSpacing[y] == exception )
        {
	  toggle = true;
        }
    }

  std::fprintf(data,"Rows,columns(#):      %8i   %8i\n",rows,columns);
  std::fprintf(data,"PixelSpacing_X,Y(mm): %8s   %8s\n",
                 pixelSpacingX,pixelSpacingY);
  std::fprintf(data,"SliceTickness(mm):    %8s\n",sliceThickness);
  std::fprintf(data,"SliceLocation(mm):    %8f\n",sliceLocation);
  std::fclose(data);

  return returnvalue;
}

G4int DicomHandler::readData(FILE *dicom,char filename2[300])
{
  G4int returnvalue = 0;
  char compressionbuf[100],maxbuf[100];
  G4int compression = 0, max = 0;
  G4int intBuffer[1000000]; 
  // intBuffer read a part of the header that is not
  //useful for this application 
  
  G4int tab[1000][1000];
  FILE* configuration = std::fopen("Data.dat","r");
  std::fscanf(configuration,"%s",compressionbuf);
  compression = atoi(compressionbuf);
  std::fscanf(configuration,"%s",maxbuf);
  max = atoi(maxbuf);
  std::fclose(configuration);

  //  READING THE PIXELS :
  G4int w = 0;
  G4int len = 0;

  if (bitAllocated == 2) // Case 16 bits :
    {
      len = rows*columns;
      for ( G4int j = 1; j <= rows; j++)
        {
	  for ( G4int i = 1; i <= columns; i++)
            {
	      w++;
	      std::fread( &intBuffer[w], 1, 2, dicom);
	      tab[j][i] = intBuffer[w];
            }
        }
    }
  else // not 16  bits :
    {
      std::printf("@@@ Error! Picture != 16 bits...\n");
      std::printf("@@@ Error! Picture != 16 bits...\n"); 
      std::printf("@@@ Error! Picture != 16 bits...\n"); 
      len = rows*columns;
      for (G4int j = 1;j <= rows;j++)
        {
	  for (G4int i = 1;i <= columns;i++)
            {
	      w++;
	      std::fread(&intBuffer[w],1,2,dicom);
	      tab[j][i] = intBuffer[w];
            }
        }
      returnvalue = 1;
    }

  // Creation of .g4 files wich contains averaged density data

  char nameProcessed[500];
  FILE* processed;

  std::sprintf(nameProcessed,"%s.g4",filename2);
  processed = std::fopen(nameProcessed,"w+");
  std::printf("### Writing of %s ###\n",nameProcessed);

  std::fprintf(processed,"%8i   %8i\n",rows,columns);
  std::fprintf(processed,"%8f   %8f\n",atof(pixelSpacingX),atof(pixelSpacingY) );
  std::fprintf(processed,"%8i\n",atoi(sliceThickness) );
  std::fprintf(processed,"%8f\n",sliceLocation);
  std::fprintf(processed,"%8i\n",compression);

  G4int compSize = 1;
  compSize = compression;
  G4int mean;
  G4bool overflow = false;
  G4int cpt=1;

  if (compSize == 1) // no compression: each pixel has a density value)
    {
      for ( G4int ww = 1; ww <= rows; ww++)
        {
	  for( G4int xx = 1; xx <= columns; xx++)
            {
	      mean = (tab[ww][xx])/1;
	      std::fprintf(processed,"%f   ",pixel2density(mean) );
            }
	  std::fprintf(processed,"\n");
        }
    }
  else
    {
      // density value is the average of a square region of
      // compression*compression pixels
      for (G4int ww = 1 ; ww <= rows ; ww = ww+compSize )
        {
	  for(G4int xx = 1 ; xx <= columns ; xx = xx+compSize )
            {
	      overflow = false;
	      G4int mean = tab[ww][xx];
	      G4int sumx = compSize-1;
	      G4int sumy = compSize-1;
	      for(  ; sumx>0; sumy--, sumx--)
                {
		  if (ww+sumy > rows|| xx+sumx > columns)
                    {
		      overflow = true;
		      break;
                    }
		  mean = mean+tab[ww+sumy][xx+sumx];
		  for (G4int m = sumx ; m>0 ; m--)
                    {
		      mean = mean+tab[ww+sumy-m][xx+sumx];
		      mean = mean+tab[ww+sumy][xx+sumx-m];
                    }
                }
	      mean = mean/(compSize*compSize);
	      cpt = 1;

	      if (overflow != true)
		std::fprintf(processed,"%f   ",pixel2density( mean) );
            }
	  std::fprintf(processed,"\n");
        }
    }
  std::fclose(processed);

  return returnvalue;
}

/*
  G4int DicomHandler::displayImage(char command[300])
  {
  //   Display DICOM images using ImageMagick
  char commandName[500];
  std::sprintf(commandName,"display  %s",command);
  std::printf(commandName);
  G4int i = system(commandName);
  return (G4int )i;
  }
*/
G4double DicomHandler::pixel2density(G4int pixel)
{
  G4double density = -1;
  G4int nbrequali = 0;
  char nbrequalibuf[100];
  G4double deltaCT = 0;
  G4double deltaDensity = 0;
  char valuedensitybuf[100][100];
  char valueCTbuf[100][100];
  G4double valuedensity[100];
  G4double valueCT[100];
  FILE* calibration;

  // CT2Density.dat contains the calibration curve to convert CT (Hounsfield) number to
  // physical density
  calibration = std::fopen("CT2Density.dat","r");
  std::fscanf(calibration,"%s",nbrequalibuf);
  nbrequali = atoi(nbrequalibuf);

  if (calibration == 0 )
    {
      std::printf("@@@ No value to transform pixels in density!\n");
      exit(1);
    }
  else // calibration != 0
    {
      for (G4int i=1;i<=nbrequali;i++) // Loop to store all the pts in CT2Density.dat
        {
	  std::fscanf(calibration,"%s %s",valueCTbuf[i-1],valuedensitybuf[i-1]);
	  valueCT[i-1] = atof(valueCTbuf[i-1]);
	  valuedensity[i-1]=atof(valuedensitybuf[i-1]);
        }
    }
  std::fclose(calibration);

  for (G4int j = 1;j<nbrequali;j++)
    {
      if ( pixel >= valueCT[j-1] && pixel < valueCT[j])
        {
	  deltaCT = valueCT[j] - valueCT[j-1];
	  deltaDensity = valuedensity[j] - valuedensity[j-1];
	  if ( pixel - valueCT[j-1] >= valueCT[j] - pixel )
            {
	      density = valuedensity[j] + ( ( valueCT[j] - pixel ) * deltaDensity/deltaCT );
            }
	  else if ( pixel - valueCT[j-1] < valueCT[j] - pixel )
            {
	      density = valuedensity[j-1] + ( ( pixel - valueCT[j-1] ) * deltaDensity/deltaCT );
            }
        }
    }

  if ( density < 0 )
    {
      std::printf("@@@ Error density = %f && Pixel = %i && deltaDensity/deltaCT = %f\n",density,pixel,deltaDensity/deltaCT);
    }

  return density;
}


void DicomHandler::checkFileFormat()
{
  std::ifstream checkData("Data.dat");
  char * oneLine = new char[101];
  G4int nbImages;

  if (!(checkData.is_open())) //Check existance of Data.dat
    {
      G4cout << "\nDicomG4 needs Data.dat :\n\tFirst line: number of image pixel for a "
	     << "voxel (G4Box)\n\tSecond line: number of images (CT slices) to "
	     << "read\n\tEach following line contains the name of a Dicom image except "
	     << "for the .dcm extension\n";
      exit(0);
    }

  checkData >> nbImages;
  checkData >> nbImages;
  G4String oneName;
  checkData.getline(oneLine,100);
  std::ifstream testExistence;
  G4bool existAlready = true;
  for (G4int rep = 0; rep < nbImages; rep++)
    { 
      checkData.getline(oneLine,100);
      oneName = oneLine;
      oneName += ".g4"; // create dicomFile.g4
      testExistence.open(oneName.data());
      if (!(testExistence.is_open()))
        {
	  existAlready = false;
	  testExistence.clear();
	  testExistence.close();
	  break;
        }
      testExistence.clear();
      testExistence.close();
    }
  checkData.close();
  delete [] oneLine;

  if ( existAlready == false ) // The files *.g4 have to be created
    {
      G4cout << "\nAll the necessary images were not found in processed form, starting "
	     << "with .dcm images\n";

      FILE* dicom;
      FILE *lecturePref;
      char compressionc[300],maxc[300];
      char name[300], inputFile[300];
      lecturePref = std::fopen("Data.dat","r");
      std::fscanf(lecturePref,"%s",compressionc);
      compression = atoi(compressionc);
      std::fscanf(lecturePref,"%s",maxc);
      max = atoi(maxc);

      for ( G4int i = 1; i <= max; i++ ) // Begin loop on filenames
        {
	  std::fscanf(lecturePref,"%s",inputFile);
	  std::sprintf(name,"%s.dcm",inputFile);
	  //  Open input file and give it to gestion_dicom :
	  std::printf("### Opening %s and reading :\n",name);
	  dicom = std::fopen(name,"rb");
	  // Reading the .dcm in two steps:
	  //      1.  reading the header
	  //	2. reading the pixel data and store the density in Moyenne.dat
	  if ( dicom != 0 )
            {
	      readHeader(dicom,inputFile);
	      readData(dicom,inputFile);
            }
	  else 
            {
	      G4cout << "\nError opening file : " << name << G4endl;
	      exit(0);
            }
	  std::fclose(dicom);
        }
      std::fclose(lecturePref);
    } 
}
