//*******************************************************
//
// DicomHandler.cc :
//	- Handling of DICM images
// 	- Reading headers and pixels
//	- Transforming pixel to density and creating *.g4
//	  files
// 	- Definitions are in DicomHandler.hh
//
// The code was written by :
//	Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
// Base on previous code by :
// 	Dragan Tubic <tdragan@gel.ulaval.ca>
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

#include "DicomHandler.hh"
#include "globals.hh"
#include "g4std/fstream"

G4int dicomHandler::readHeader(FILE *dicom, char filename2[300])
{
  G4int returnvalue=0;
  // The first 128 bytes are not important
  fread( buffer, 1, 128, dicom );

  // Reads the "DICM" letters
  fread( buffer, 1, 4, dicom );

  // Read information up to the pixel data
  // note: it should be a while instead of a for
  for (G4int i=0;i<=100000000;i++)
    {
      //Reading groups and elements :
      fread(&read_group_id,1,2,dicom);
      fread(&read_element_id,1,2,dicom);

      if (read_group_id == 0x7FE0) // beginning of the pixels
        {
	  fread( buffer, 1, 2,dicom); // Skip 2 reserved bytes
	  break;
        }

      fread(&element_length,1,2,dicom);

      // If value representation (VR) is OB, OW, SQ, UN, the next length is 32 bits
      if (element_length == 16975 || element_length == 22351 || element_length == 20819 || element_length == 20053)
        {
	  //skip 2 reserved "bytes"
	  fread( buffer, 1, 2,dicom); // Skip 2 reserved bytes
	  fread(&element_length3, 4, 1, dicom);// Reading length of the information
	  fread(&value[i],element_length3,1,dicom);// Reading the information with
	  // (BIG) buffer : "value"
	  // Creating a tag to be identified afterward
	  tag_dictionnary=read_group_id*0x10000 + read_element_id;
        }
      else  // lenght is 16 bits :
        {
	  fread(&element_length2,1,2,dicom);
	  fread(&value[i],element_length2,1,dicom);
	  tag_dictionnary=read_group_id*0x10000 + read_element_id;
        }

      if (tag_dictionnary == 0x00280010 ) // Number of Rows
        {
	  rows=*(G4int*)&value[i];
	  printf("[0x00280010] Rows -> %i\n",rows);
        }
      if (tag_dictionnary == 0x00280011 ) // Number of columns
        {
	  columns=*(G4int*)&value[i];
	  printf("[0x00280011] Columns -> %i\n",columns);
        }
      if (tag_dictionnary == 0x00280102 ) // High bits  ( not used )
        {
	  high_bits=*(G4int*)&value[i];
	  printf("[0x00280102] High bits -> %i\n",high_bits);
        }
      if (tag_dictionnary == 0x00280100 )  // Bits allocated ( not used )
        {
	  bits_allocated=*(G4int*)&value[i];
	  printf("[0x00280100] Bits allocated -> %i\n",bits_allocated);
	  bits_allocated=(bits_allocated)/8;
        }
      if (tag_dictionnary == 0x00280101 )  //  Bits stored ( not used )
        {
	  bits_stored=*(G4int*)&value[i];
	  printf("[0x00280101] Bits stord -> %i\n",bits_stored);
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00281053)  //  Rescale slope ( not used )
        {
	  rescale_slope= atoi( value[i] );
	  printf("[0x00281053] Rescale Slope -> %i\n",  atoi( value[i] ) );
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00281052 )  // Rescalse intercept ( not used )
        {
	  rescale_intercept = atoi( value[i] );
	  printf("[0x00281052] Rescale Intercept -> %i\n",  atoi( value[i] ) );
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00280103 )
        {
	  //  Pixel representation ( functions not design to read signed bits )
	  printf("[0x00280103] Pixel Representation -> %i\n",  atoi( value[i] ) );
	  if ( atoi(value[i]) == 1 )
            {
	      printf("### PIXEL REPRESENTATION = 1, BITS ARE SIGNED, ");
	      printf("DICOM READING SCAN FOR UNSIGNED VALUE, POSSIBLE ");
	      printf("ERROR !!!!!! -> \n");
            }
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00080008 ) //  Image type ( not used )
        {
	  printf("[0x00080008] Image Types -> %s\n",value[i]);
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00283000 )   //  Modality LUT Sequence ( not used )
        {
	  printf("[0x00283000] Modality LUT Sequence SQ 1 -> %s\n",value[i]);
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00283002 )  // LUT Descriptor ( not used )
        {
	  printf("[0x00283002] LUT Descriptor US or SS 3 -> %s\n",value[i]);
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00283003 )  // LUT Explanation ( not used )
        {
	  printf("[0x00283003] LUT Explanation LO 1 -> %s\n",value[i]);
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00283004 )  // Modality LUT ( not used )
        {
	  printf("[0x00283004] Modality LUT Type LO 1 -> %s\n",value[i]);
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00283006 )  // LUT Data ( not used )
        {
	  printf("[0x00283006] LUT Data US or SS -> %s\n",value[i]);
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00283010 )  // VOI LUT ( not used )
        {
	  printf("[0x00283010] VOI LUT Sequence SQ 1 -> %s\n",value[i]);
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00280120 )  // Pixel Padding Value ( not used )
        {
	  printf("[0x00280120] Pixel Padding Value US or SS 1 -> %s\n",value[i]);
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00280030 )  // Pixel Spacing
        {
	  printf("[0x00280030] Pixel Spacing (mm) -> %s\n",value[i]);
	  sprintf(pixel_spacing,"%s",value[i]);//pixel_spacing=value[i];
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00200037 )  // Image Orientation ( not used )
        {
	  printf("[0x00200037] Image Orientation (Patient) -> %s\n",value[i]);
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00200032 )  // Image Position ( not used )
        {
	  printf("[0x00200032] Image Position (Patient,mm) -> %s\n",value[i]);
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00180050 )  // Slice Tickness
        {
	  printf("[0x00180050] Slice Tickness (mm) -> %s\n",value[i]);
	  sprintf(slice_tickness,"%s",value[i]);//slice_tickness=value[i];
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00201041 )  // Slice Location
        {
	  printf("[0x00201041] Slice Location -> %s\n",value[i]);
	  slice_location=atof(value[i]);
	  bits_stored=(bits_stored)/8;
        }
      if (tag_dictionnary == 0x00280004 )  // Photometric Interpretation ( not used )
        {
	  printf("[0x00280004] Photometric Interpretation -> %s\n",value[i]);
	  bits_stored=(bits_stored)/8;
        }


    }

  // Creating files to store information
  char compressionbuf[100],maxbuf[100];
  char filename[300];
  G4int compression=0;
  G4int  max=0;
  FILE* configuration;

  configuration=fopen("Data.dat","r");
  if ( configuration != NULL )
    {
      fscanf(configuration,"%s",compressionbuf);
      compression=atoi(compressionbuf);
      fscanf(configuration,"%s",maxbuf);
      max=atoi(maxbuf);
      fclose(configuration);
    }
  else
    {
      printf("### WARNING, file Data.dat not here !!!\n");
      exit(1);
    }

  sprintf(filename,"%s.dat",filename2);
  data = fopen(filename,"w+");
  // Note: the .dat files contain basic information on the images.

  char exception = '\\';
  G4bool toggle=false;
  G4int z=0;
  for (G4int y=0;y<=300;y++)
    {
      if ( pixel_spacing[y] != exception )
        {
	  if (toggle == false)
	    pixel_spacing_X[y] = pixel_spacing[y];
	  if (toggle == true)
            {
	      pixel_spacing_Y[z] = pixel_spacing[y];
	      z++;
            }
        }
      else if ( pixel_spacing[y] == exception )
        {
	  toggle=true;
        }
    }

  fprintf(data,"Rows,columns(#):      %8i   %8i\n",rows,columns);
  fprintf(data,"PixelSpacing_X,Y(mm): %8s   %8s\n",pixel_spacing_X,pixel_spacing_Y);
  fprintf(data,"SliceTickness(mm):    %8s\n",slice_tickness);
  fprintf(data,"SliceLocation(mm):    %8f\n",slice_location);
  fclose(data);

  return returnvalue;
}

G4int dicomHandler::readData(FILE *dicom,	char filename2[300])
{
  G4int returnvalue=0;
  char compressionbuf[100],maxbuf[100];
  G4int compression=0, max=0;

  FILE* configuration=fopen("Data.dat","r");
  fscanf(configuration,"%s",compressionbuf);
  compression=atoi(compressionbuf);
  fscanf(configuration,"%s",maxbuf);
  max=atoi(maxbuf);
  fclose(configuration);

  //  READING THE PIXELS :
  G4int w=0;
  if (bits_allocated == 2) // Case 16 bits :
    {
      len = rows*columns;
      for (G4int j=1;j<=rows;j++)
        {
	  for (G4int i=1;i<=columns;i++)
            {
	      w++;
	      fread(&Int_Buffer[w],1,2,dicom);
	      tab[j][i]=Int_Buffer[w];
            }
        }
    }
  else // not 16  bits :
    {
      printf("@@@ Error! Picture != 16 bits...\n"); // Expect the program to CRASH !!!
      printf("@@@ Error! Picture != 16 bits...\n"); // Expect the program to CRASH !!!
      printf("@@@ Error! Picture != 16 bits...\n"); // Expect the program to CRASH !!!
      // We still try to read the image as if it were 16 bits
      len = rows*columns;
      for (G4int j=1;j<=rows;j++)
        {
	  for (G4int i=1;i<=columns;i++)
            {
	      w++;
	      fread(&Int_Buffer[w],1,2,dicom);
	      tab[j][i]=Int_Buffer[w];
            }
        }
      returnvalue=1;
    }

  // Creation of .g4 files wich contains averaged density data

  char nameProcessed[500];
  FILE* processed;

  sprintf(nameProcessed,"%s.g4",filename2);
  processed = fopen(nameProcessed,"w+");
  printf("### Writing of %s ###\n",nameProcessed);

  fprintf(processed,"%8i   %8i\n",rows,columns);
  fprintf(processed,"%8f   %8f\n",atof(pixel_spacing_X),atof(pixel_spacing_Y) );
  fprintf(processed,"%8i\n",atoi(slice_tickness) );
  fprintf(processed,"%8f\n",slice_location);
  fprintf(processed,"%8i\n",compression);

  G4int compSize = 1;
  compSize=compression;
  G4int mean;
  G4bool overflow=false;
  G4int cpt=1;

  if (compSize == 1) // no compression: each pixel has a density value)
    {
      for (G4int ww =1;ww<=rows; ww++)
        {
	  for(G4int xx =1 ;xx<=columns;xx++)
            {
	      mean=(tab[ww][xx])/1;
	      fprintf(processed,"%f   ",pixel2density(mean) );
            }
	  fprintf(processed,"\n");
        }
    }
  else
    {
      // density value is the average of a square region of
      // compression*compression pixels
      for (G4int ww =1 ; ww<=rows ; ww=ww+compSize )
        {
	  for(G4int xx =1 ; xx<=columns ; xx=xx+compSize )
            {
	      overflow=false;
	      G4int mean=tab[ww][xx];
	      G4int sumx=compSize-1;
	      G4int sumy=compSize-1;
	      for(  ; sumx>0; sumy--, sumx--)
                {
		  if (ww+sumy > rows|| xx+sumx > columns)
                    {
		      overflow=true;
		      break;
                    }
		  mean=mean+tab[ww+sumy][xx+sumx];
		  for (G4int m=sumx ; m>0 ; m--)
                    {
		      mean=mean+tab[ww+sumy-m][xx+sumx];
		      mean=mean+tab[ww+sumy][xx+sumx-m];
                    }
                }
	      mean=mean/(compSize*compSize);
	      cpt=1;

	      if (overflow != true)
		fprintf(processed,"%f   ",pixel2density( mean) );
            }
	  fprintf(processed,"\n");
        }
    }
  fclose(processed);


  return returnvalue;
}

G4int dicomHandler::displayImage(char command[300])
{
  //   Display DICOM images using ImageMagick
  char commandName[500];
  sprintf(commandName,"display  %s",command);
  printf(commandName);
  G4int i=system(commandName);
  return (G4int )i;
}

G4double dicomHandler::pixel2density(G4int pixel)
{
  G4double density=-1;
  G4int nbrequali=0;
  char nbrequalibuf[100];
  G4double deltaCT=0;
  G4double deltaDensity=0;
  char valuedensitybuf[100][100];
  char valueCTbuf[100][100];
  G4double valuedensity[100];
  G4double valueCT[100];
  FILE* calibration;

  // CT2Density.dat contains the calibration curve to convert CT (Hounsfield) number to
  // physical density
  calibration=fopen("CT2Density.dat","r");
  fscanf(calibration,"%s",nbrequalibuf);
  nbrequali=atoi(nbrequalibuf);

  if (calibration == NULL )
    {
      printf("@@@ No value to transform pixels in density!\n");
      exit(1);
    }
  else // calibration != NULL
    {
      for (G4int i=1;i<=nbrequali;i++) // Loop to store all the pts in CT2Density.dat
        {
	  fscanf(calibration,"%s %s",valueCTbuf[i-1],valuedensitybuf[i-1]);
	  valueCT[i-1]=atof(valueCTbuf[i-1]);
	  valuedensity[i-1]=atof(valuedensitybuf[i-1]);
        }
    }
  fclose(calibration);

  for (G4int j=1;j<nbrequali;j++)
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
      printf("@@@ Error density = %f && Pixel = %i && deltaDensity/deltaCT = %f\n",density,pixel,deltaDensity/deltaCT);
    }

  return density;
}


void dicomHandler::checkFileFormat()
{
  G4cout << "\n\n\n#+#+#+#+#+#+#+#+#+#+#+#+#+##+#+#+#+#+#+#+#+#+#\n"
	 << "  Checking for previously processed image set ...";

  G4std::ifstream checkData("Data.dat");
  char * oneLine = new char[101];
  G4int nbImages;

  if (!(checkData.is_open()))
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
  G4std::ifstream testExistence;
  G4bool existAlready = true;
  for (G4int rep=0; rep<nbImages; rep++)
    {
      checkData.getline(oneLine,100);
      oneName = oneLine;
      oneName += ".g4";
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
      FILE *lecturepref;
      char compressionc[300],maxc[300];
      lecturepref=fopen("Data.dat","r");
      fscanf(lecturepref,"%s",compressionc);
      compression=atoi(compressionc);
      fscanf(lecturepref,"%s",maxc);
      max=atoi(maxc);

      for (G4int i=1;i<=max;i++) // Begin loop on filenames
        {
	  fscanf(lecturepref,"%s",name_in_file);
	  sprintf(name,"%s.dcm",name_in_file);
	  //  Open input file and give it to gestion_dicom :
	  printf("### Opening %s and reading :\n",name);
	  dicom = fopen(name,"rb");
	  // Reading the .dcm in two steps:
	  //      1.  reading the header
	  //	2. reading the pixel data and store the density in Moyenne.dat
	  if ( dicom != NULL )
            {
	      readHeader(dicom,name_in_file);
	      readData(dicom,name_in_file);
            }
	  else //if ( dicom == NULL )
            {
	      G4cout << "\nError opening file : " << name << G4endl;
	      exit(0);
            }
	  fclose(dicom);
        }// End loop on filenames
      fclose(lecturepref);
    } // files .g4 built
  else
    G4cout << "\nProcessed images found, a bit of time saved\n";

  G4cout << "#+#+#+#+#+#+#+#+#+#+#+#+#+##+#+#+#+#+#+#+#+#+#\n\n\n";
}
