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
// + Universit.A�� Laval, Qu��bec (QC) Canada
//*******************************************************
//
//*******************************************************
//
//*******************************************************
//
// DicomHandler.cc :
//	- Handling of DICM images
// 	- Reading headers and pixels
//	- Transforming pixel to density and creating *.g4dcm
//	  files
// 	- Definitions are in DicomHandler.hh
//*******************************************************
//
// The code is part of the DICOM extended example and it was modified by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#include "DicomHandler.hh"
#include "globals.hh"
#include "G4ios.hh"
#include <fstream>

#include <cctype>
#include <cstring>


DicomHandler::DicomHandler()
    : DATABUFFSIZE(8192), LINEBUFFSIZE(128), FILENAMESIZE(512),
      compression(0), nFiles(0), rows(0), columns(0),
      bitAllocated(0), maxPixelValue(0), minPixelValue(0),
      pixelSpacingX(0.), pixelSpacingY(0.),
      sliceThickness(0.), sliceLocation(0.),
      rescaleIntercept(0), rescaleSlope(0),
      littleEndian(true), implicitEndian(false),
      pixelRepresentation(0) {
    ;
}

DicomHandler::~DicomHandler() {
    ;
}

G4int DicomHandler::ReadFile(FILE *dicom, char * filename2)
{
  G4cout << " ReadFile " << filename2 << G4endl;
    G4int returnvalue = 0;
    char * buffer = new char[LINEBUFFSIZE];

    implicitEndian = false;
    littleEndian = true;

    std::fread( buffer, 1, 128, dicom ); // The first 128 bytes 
                                         //are not important
    // Reads the "DICOM" letters
    std::fread( buffer, 1, 4, dicom );
    // if there is no preamble, the FILE pointer is rewinded.
    if(std::strncmp("DICM", buffer, 4) != 0) {
        std::fseek(dicom, 0, SEEK_SET);
        implicitEndian = true;
    }

    short readGroupId;    // identify the kind of input data 
    short readElementId;  // identify a particular type information
    short elementLength2; // deal with element length in 2 bytes
    //unsigned int elementLength4; // deal with element length in 4 bytes
    G4int elementLength4; // deal with element length in 4 bytes

    char * data = new char[DATABUFFSIZE];

 
    // Read information up to the pixel data
    while(true) {

	//Reading groups and elements :
	readGroupId = 0;
	readElementId = 0;
	// group ID
	std::fread(buffer, 2, 1, dicom);
	GetValue(buffer, readGroupId);
	// element ID
	std::fread(buffer, 2, 1, dicom);
	GetValue(buffer, readElementId);

	// Creating a tag to be identified afterward
	G4int tagDictionary = readGroupId*0x10000 + readElementId;

      // VR or element length
	std::fread(buffer,2,1,dicom);
	GetValue(buffer, elementLength2);
	  
 	// If value representation (VR) is OB, OW, SQ, UN, added OF and UT 
	//the next length is 32 bits
	if((elementLength2 == 0x424f ||  // "OB"
	    elementLength2 == 0x574f ||  // "OW"
	    elementLength2 == 0x464f ||  // "OF"
	    elementLength2 == 0x5455 ||  // "UT"
	    elementLength2 == 0x5153 || //  "SQ"
	    elementLength2 == 0x4e55) && // "UN"
	   !implicitEndian ) {           // explicit VR

	    std::fread(buffer, 2, 1, dicom); // Skip 2 reserved bytes

	    // element length
	    std::fread(buffer, 4, 1, dicom);
	    GetValue(buffer, elementLength4);

	    // beginning of the pixels
	    if(tagDictionary == 0x7FE00010) break;
	    
	    if(elementLength4 ==(G4int) 0xFFFFFFFF) G4cerr << "Too difficult for me!" << G4endl; // VR = SQ with undefined length
	    
	    // Reading the information with data
	    std::fread(data, elementLength4,1,dicom);

		
	}  else { 

		if(!implicitEndian || readGroupId == 2) {  //  explicit with VR different than previous ones
	    	
		  //G4cout << "Reading  DICOM files with Explicit VR"<< G4endl;
		  // element length (2 bytes)
		  std::fread(buffer, 2, 1, dicom);
		  GetValue(buffer, elementLength2);
		  elementLength4 = elementLength2;
		  
		  if(tagDictionary == 0x7FE00010) break; // beginning of the pixels
	          
		  std::fread(data, elementLength4, 1, dicom);
	        
		} else { 				 // Implicit VR

                  //G4cout << "Reading  DICOM files with Implicit VR"<< G4endl;
   
		  // element length (4 bytes)
		  if(std::fseek(dicom, -2, SEEK_CUR) != 0) {
		      G4cerr << "[DicomHandler] fseek failed" << G4endl;
		      exit(-10);}

		  std::fread(buffer, 4, 1, dicom);
		  GetValue(buffer, elementLength4);

	          //G4cout <<  std::hex<< elementLength4 << G4endl;
	      
	          if(elementLength4 ==(G4int) 0xFFFFFFFF) {
	          short momreadGroupId;
		  short momreadElementId;
		  unsigned int momelementLength4 ; // deal with element length in 4 bytes
      		  G4cout << "Trying to read nested items!" << G4endl; // VR = SQ with undefined length
		  while(true){
			std::fread(buffer, 2, 1, dicom);
			GetValue(buffer, momreadGroupId);
			std::fread(buffer, 2, 1, dicom);
			GetValue(buffer, momreadElementId);
			//G4int momtagDictionary = momreadGroupId*0x10000 + momreadElementId;
			//G4cout << "TAG"<< std::hex << momtagDictionary << G4endl;
			std::fread(buffer, 4, 1, dicom);
			GetValue(buffer, momelementLength4);
			//G4cout << "Length"<< std::hex << momelementLength4 << G4endl;
			if(momelementLength4 == 0x00000000 || momelementLength4 == 0xFFFFFFFF ) break;
                	std::fread(buffer, momelementLength4, 1, dicom);
			}
		  } else {
	          if(tagDictionary == 0x7FE00010) break; // beginning of the pixels
	          
		  std::fread(data, elementLength4, 1, dicom);
	        
		 } 
		      
	       } 
	}

	// NULL termination
	data[elementLength4] = '\0';

	// analyzing information 
	GetInformation(tagDictionary, data);
    }

    // Creating files to store information
    std::ofstream foutG4DCM;
	char * fullFileName=new char[1024];
        strcpy(fullFileName,this->dicomDirectory);
        strcat(fullFileName, filename2);
    G4String fnameG4DCM = G4String(fullFileName) + ".g4dcm";
    foutG4DCM.open(fnameG4DCM);
    G4cout << "### Writing of " << fnameG4DCM << " ### " << G4endl;

    foutG4DCM << fMaterialIndices.size() << G4endl;
    //--- Write materials
    unsigned int ii = 0;
    std::map<G4float,G4String>::const_iterator ite;
    for( ite = fMaterialIndices.begin(); ite != fMaterialIndices.end(); ite++, ii++ ){
      foutG4DCM << ii << " " << (*ite).second << G4endl;
    }
    //--- Write number of voxels (assume only one voxel in Z)
    foutG4DCM << rows/compression << " " << columns/compression << " 1 " << G4endl;
    //--- Write minimum and maximum extensions
    foutG4DCM << -pixelSpacingX*rows/2 << " " << pixelSpacingX*rows/2 << G4endl;
    foutG4DCM << -pixelSpacingY*columns/2 << " " << pixelSpacingY*columns/2 << G4endl;
    foutG4DCM << sliceLocation-sliceThickness/2. << " " << sliceLocation+sliceThickness/2. << G4endl;
    //    foutG4DCM << compression << G4endl;
    
    ReadData( dicom, fullFileName );
    
    StoreData( foutG4DCM );

    foutG4DCM.close();

    //
    delete [] buffer;
    delete [] data;

    return returnvalue;
}

//
void DicomHandler::GetInformation(G4int & tagDictionary, char * data) {
    if(tagDictionary == 0x00280010 ) { // Number of Rows
	GetValue(data, rows);
	std::printf("[0x00280010] Rows -> %i\n",rows);

    } else if(tagDictionary == 0x00280011 ) { // Number of columns
	GetValue(data, columns);
	std::printf("[0x00280011] Columns -> %i\n",columns);

    } else if(tagDictionary == 0x00280102 ) { // High bits  ( not used )
	short highBits;
	GetValue(data, highBits);
	std::printf("[0x00280102] High bits -> %i\n",highBits);

    } else if(tagDictionary == 0x00280100 ) { // Bits allocated
	GetValue(data, bitAllocated);
	std::printf("[0x00280100] Bits allocated -> %i\n", bitAllocated);

    } else if(tagDictionary == 0x00280101 ) { //  Bits stored ( not used )
	short bitStored;
	GetValue(data, bitStored);
	std::printf("[0x00280101] Bits stored -> %i\n",bitStored);

    } else if(tagDictionary == 0x00280106 ) { //  Min. pixel value
	GetValue(data, minPixelValue);
	std::printf("[0x00280106] Min. pixel value -> %i\n", minPixelValue);

    } else if(tagDictionary == 0x00280107 ) { //  Max. pixel value
	GetValue(data, maxPixelValue);
	std::printf("[0x00280107] Max. pixel value -> %i\n", maxPixelValue);

    } else if(tagDictionary == 0x00281053) { //  Rescale slope
	rescaleSlope = atoi(data);
	std::printf("[0x00281053] Rescale Slope -> %d\n", rescaleSlope);

    } else if(tagDictionary == 0x00281052 ) { // Rescalse intercept
	rescaleIntercept = atoi(data);
	std::printf("[0x00281052] Rescale Intercept -> %d\n", rescaleIntercept );

    } else if(tagDictionary == 0x00280103 ) {
	//  Pixel representation ( functions not design to read signed bits )
	pixelRepresentation = atoi(data); // 0: unsigned  1: signed 
	std::printf("[0x00280103] Pixel Representation -> %i\n", pixelRepresentation);
	if(pixelRepresentation == 1 ) {
	    std::printf("### PIXEL REPRESENTATION = 1, BITS ARE SIGNED, ");
	    std::printf("DICOM READING SCAN FOR UNSIGNED VALUE, POSSIBLE ");
	    std::printf("ERROR !!!!!! -> \n");
	}

    } else if(tagDictionary == 0x00080006 ) { //  Modality
	std::printf("[0x00080006] Modality -> %s\n", data);

    } else if(tagDictionary == 0x00080070 ) { //  Manufacturer
	std::printf("[0x00080070] Manufacturer -> %s\n", data);

    } else if(tagDictionary == 0x00080080 ) { //  Institution Name
	std::printf("[0x00080080] Institution Name -> %s\n", data);

    } else if(tagDictionary == 0x00080081 ) { //  Institution Address
	std::printf("[0x00080081] Institution Address -> %s\n", data);

    } else if(tagDictionary == 0x00081040 ) { //  Institution Department Name
	std::printf("[0x00081040] Institution Department Name -> %s\n", data);

    } else if(tagDictionary == 0x00081090 ) { //  Manufacturer's Model Name
	std::printf("[0x00081090] Manufacturer's Model Name -> %s\n", data);

    } else if(tagDictionary == 0x00181000 ) { //  Device Serial Number
	std::printf("[0x00181000] Device Serial Number -> %s\n", data);

    } else if(tagDictionary == 0x00080008 ) { //  Image type ( not used )
	std::printf("[0x00080008] Image Types -> %s\n", data);
	    
    } else if(tagDictionary == 0x00283000 ) { //  Modality LUT Sequence ( not used )
	std::printf("[0x00283000] Modality LUT Sequence SQ 1 -> %s\n", data);

    } else if(tagDictionary == 0x00283002 ) { // LUT Descriptor ( not used )
	std::printf("[0x00283002] LUT Descriptor US or SS 3 -> %s\n", data);

    } else if(tagDictionary == 0x00283003 ) { // LUT Explanation ( not used )
	std::printf("[0x00283003] LUT Explanation LO 1 -> %s\n", data);

    } else if(tagDictionary == 0x00283004 ) { // Modality LUT ( not used )
	std::printf("[0x00283004] Modality LUT Type LO 1 -> %s\n", data);

    } else if(tagDictionary == 0x00283006 ) { // LUT Data ( not used )
	std::printf("[0x00283006] LUT Data US or SS -> %s\n", data);

    } else if(tagDictionary == 0x00283010 ) { // VOI LUT ( not used )
	std::printf("[0x00283010] VOI LUT Sequence SQ 1 -> %s\n", data);

    } else if(tagDictionary == 0x00280120 ) { // Pixel Padding Value ( not used )
	std::printf("[0x00280120] Pixel Padding Value US or SS 1 -> %s\n", data);

    } else if(tagDictionary == 0x00280030 ) { // Pixel Spacing
      G4String datas(data);
      int iss = datas.find('\\');
      pixelSpacingX = atof( datas.substr(0,iss).c_str() );
      pixelSpacingY = atof( datas.substr(iss+2,datas.length()).c_str() );

    } else if(tagDictionary == 0x00200037 ) { // Image Orientation ( not used )
	std::printf("[0x00200037] Image Orientation (Patient) -> %s\n", data);

    } else if(tagDictionary == 0x00200032 ) { // Image Position ( not used )
	std::printf("[0x00200032] Image Position (Patient,mm) -> %s\n", data);

    } else if(tagDictionary == 0x00180050 ) { // Slice Thickness
	sliceThickness = atof(data);
	std::printf("[0x00180050] Slice Thickness (mm) -> %f\n", sliceThickness);

    } else if(tagDictionary == 0x00201041 ) { // Slice Location
	sliceLocation = atof(data);
	std::printf("[0x00201041] Slice Location -> %f\n", sliceLocation);

    } else if(tagDictionary == 0x00280004 ) { // Photometric Interpretation ( not used )
	std::printf("[0x00280004] Photometric Interpretation -> %s\n", data);

    } else if(tagDictionary == 0x00020010) { // Endian
	if(strcmp(data, "1.2.840.10008.1.2") == 0)
	    implicitEndian = true;
	else if(strncmp(data, "1.2.840.10008.1.2.2", 19) == 0)
	    littleEndian = false;
	//else 1.2.840..10008.1.2.1 (explicit little endian)
		   
	std::printf("[0x00020010] Endian -> %s\n", data);
    }

    // others
    else {
	std::printf("[0x%x] -> %s\n", tagDictionary, data);

    }

}

void DicomHandler::StoreData(std::ofstream& foutG4DCM) 
{
  G4int mean=0;
  G4double density;
  G4bool overflow = false;
  G4int cpt=1;


  //----- Print indices of material 
  if(compression == 1) { // no compression: each pixel has a density value)
    for( G4int ww = 0; ww < rows; ww++) {
      for( G4int xx = 0; xx < columns; xx++) {
	mean = tab[ww][xx];
	density = Pixel2density(mean);
	foutG4DCM << GetMaterialIndex( density ) << " "; 
      }
      foutG4DCM << G4endl;
    }
    
  } else {
    // density value is the average of a square region of
    // compression*compression pixels
    for(G4int ww = 0; ww < rows ;ww += compression ) {
      for(G4int xx = 0; xx < columns ;xx +=compression ) {
	overflow = false;
	mean = 0;
	for(int sumx = 0; sumx < compression; sumx++) {
	  for(int sumy = 0; sumy < compression; sumy++) {
	    if(ww+sumy >= rows || xx+sumx >= columns) overflow = true;
	    mean += tab[ww+sumy][xx+sumx];
	  }
	  if(overflow) break;
	}
	mean /= compression*compression;
	cpt = 1;
	
	if(!overflow) {
	  G4double density = Pixel2density(mean);
	  foutG4DCM << GetMaterialIndex( density ) << " ";
	}
      }
      foutG4DCM << G4endl;
    }

  }

  //----- Print densities
  if(compression == 1) { // no compression: each pixel has a density value)
    for( G4int ww = 0; ww < rows; ww++) {
      for( G4int xx = 0; xx < columns; xx++) {
	mean = tab[ww][xx];
	density = Pixel2density(mean);
	foutG4DCM << density << " ";
	if( xx%8 == 3 ) foutG4DCM << G4endl; // just for nicer reading
      }
    }
    
  } else {
    // density value is the average of a square region of
    // compression*compression pixels
    for(G4int ww = 0; ww < rows ;ww += compression ) {
      for(G4int xx = 0; xx < columns ;xx +=compression ) {
	overflow = false;
	mean = 0;
	for(int sumx = 0; sumx < compression; sumx++) {
	  for(int sumy = 0; sumy < compression; sumy++) {
	    if(ww+sumy >= rows || xx+sumx >= columns) overflow = true;
	    mean += tab[ww+sumy][xx+sumx];
	  }
	  if(overflow) break;
	}
	mean /= compression*compression;
	cpt = 1;
	
	if(!overflow) {
	  G4double density = Pixel2density(mean);
	  foutG4DCM << density  << " ";
	  if( xx/compression%8 == 3 ) foutG4DCM << G4endl; // just for nicer reading
	}
      }
    }

  }

}

void DicomHandler::ReadMaterialIndices( std::ifstream& finData)
{
  unsigned int nMate;
  G4String mateName;
  G4float densityMax;
  finData >> nMate;
  if( finData.eof() ) return;

  G4cout << " ReadMaterialIndices " << nMate << G4endl;
  for( unsigned int ii = 0; ii < nMate; ii++ ){
    finData >> mateName >> densityMax;
    fMaterialIndices[densityMax] = mateName;
    G4cout << ii << " ReadMaterialIndices " << mateName << " " << densityMax << G4endl;
  }

}

unsigned int DicomHandler::GetMaterialIndex( G4float density )
{
 std::map<G4float,G4String>::reverse_iterator ite;
  G4int ii = fMaterialIndices.size();
  for( ite = fMaterialIndices.rbegin(); ite != fMaterialIndices.rend(); ite++, ii-- ) {
    if( density >= (*ite).first ) { 
      break;
    }
  }
//     G4cout << " GetMaterialIndex " << density << " = " << ii << G4endl;
//     G4cout << density << " " << ii << G4endl;
  return  ii;

}

//
G4int DicomHandler::ReadData(FILE *dicom,char * filename2)
{
    G4int returnvalue = 0;

    //  READING THE PIXELS :
    G4int w = 0;
    G4int len = 0;
    
    tab = new G4int*[rows];
    for ( G4int i = 0; i < rows; i ++ ) {
      tab[i] = new G4int[columns];
    }

    if(bitAllocated == 8) { // Case 8 bits :

	std::printf("@@@ Error! Picture != 16 bits...\n");
	std::printf("@@@ Error! Picture != 16 bits...\n"); 
	std::printf("@@@ Error! Picture != 16 bits...\n"); 

	unsigned char ch = 0;

	len = rows*columns;
	for(G4int j = 0; j < rows; j++) {
	    for(G4int i = 0; i < columns; i++) {
		w++;
		std::fread( &ch, 1, 1, dicom);
		tab[j][i] = ch*rescaleSlope + rescaleIntercept;
	    }
	}
	returnvalue = 1;

    } else { //  from 12 to 16 bits :
	char sbuff[2];
	short pixel;
	len = rows*columns;
	for( G4int j = 0; j < rows; j++) {
	    for( G4int i = 0; i < columns; i++) {
		w++;
		std::fread(sbuff, 2, 1, dicom);
		GetValue(sbuff, pixel);
		tab[j][i] = pixel*rescaleSlope + rescaleIntercept;
	    }
	}
    }

    // Creation of .g4 files wich contains averaged density data
    char * nameProcessed = new char[FILENAMESIZE];
    FILE* fileOut;

    std::sprintf(nameProcessed,"%s.g4dcmb",filename2);
    fileOut = std::fopen(nameProcessed,"w+b");
    std::printf("### Writing of %s ###\n",nameProcessed);

    unsigned int nMate = fMaterialIndices.size();
    std::fwrite(&nMate, sizeof(unsigned int), 1, fileOut);
    //--- Write materials
    std::map<G4float,G4String>::const_iterator ite;
    for( ite = fMaterialIndices.begin(); ite != fMaterialIndices.end(); ite++ ){
      G4String mateName = (*ite).second;
      for( G4int ii = (*ite).second.length(); ii < 40; ii++ ) {
	mateName += " ";
      }	 //mateName = const_cast<char*>(((*ite).second).c_str());

      const char* mateNameC = mateName.c_str();
      std::fwrite(mateNameC, sizeof(char),40, fileOut);
    }

    unsigned int rowsC = rows/compression;
    unsigned int columnsC = columns/compression;
    unsigned int planesC = 1;
    G4float pixelLocationXM = -pixelSpacingX*rows/2.;
    G4float pixelLocationXP = pixelSpacingX*rows/2.;
    G4float pixelLocationYM = -pixelSpacingY*rows/2.;
    G4float pixelLocationYP = pixelSpacingY*rows/2.;
    G4float sliceLocationZM = sliceLocation-sliceThickness/2.;
    G4float sliceLocationZP = sliceLocation+sliceThickness/2.;
    //--- Write number of voxels (assume only one voxel in Z)
    std::fwrite(&rowsC, sizeof(unsigned int), 1, fileOut);
    std::fwrite(&columnsC, sizeof(unsigned int), 1, fileOut);
    std::fwrite(&planesC, sizeof(unsigned int), 1, fileOut);
    //--- Write minimum and maximum extensions
    std::fwrite(&pixelLocationXM, sizeof(G4float), 1, fileOut);
    std::fwrite(&pixelLocationXP, sizeof(G4float), 1, fileOut);
    std::fwrite(&pixelLocationYM, sizeof(G4float), 1, fileOut);
    std::fwrite(&pixelLocationYP, sizeof(G4float), 1, fileOut);
    std::fwrite(&sliceLocationZM, sizeof(G4float), 1, fileOut);
    std::fwrite(&sliceLocationZP, sizeof(G4float), 1, fileOut);
    //    std::fwrite(&compression, sizeof(unsigned int), 1, fileOut);

    std::printf("%8i   %8i\n",rows,columns);
    std::printf("%8f   %8f\n",pixelSpacingX,pixelSpacingY);
    std::printf("%8f\n", sliceThickness);
    std::printf("%8f\n", sliceLocation);
    std::printf("%8i\n", compression);

    G4int compSize = compression;
    G4int mean;
    G4float density;
    G4bool overflow = false;
    G4int cpt=1;

    //----- Write index of material for each pixel
    if(compSize == 1) { // no compression: each pixel has a density value)
      for( G4int ww = 0; ww < rows; ww++) {
	for( G4int xx = 0; xx < columns; xx++) {
	  mean = tab[ww][xx];
	  density = Pixel2density(mean);
	  unsigned int mateID = GetMaterialIndex( density );
	  std::fwrite(&mateID, sizeof(unsigned int), 1, fileOut);
	}
      }

    } else {
      // density value is the average of a square region of
      // compression*compression pixels
      for(G4int ww = 0; ww < rows ;ww += compSize ) {
	for(G4int xx = 0; xx < columns ;xx +=compSize ) {
	  overflow = false;
	  mean = 0;
	  for(int sumx = 0; sumx < compSize; sumx++) {
	    for(int sumy = 0; sumy < compSize; sumy++) {
	      if(ww+sumy >= rows || xx+sumx >= columns) overflow = true;
	      mean += tab[ww+sumy][xx+sumx];
	    }
	    if(overflow) break;
	  }
	  mean /= compSize*compSize;
	  cpt = 1;
	  
	  if(!overflow) {
	    density = Pixel2density(mean);
	    unsigned int mateID = GetMaterialIndex( density );
	    std::fwrite(&mateID, sizeof(unsigned int), 1, fileOut);
	  }
	}
	
      }
    }

    //----- Write density for each pixel
    if(compSize == 1) { // no compression: each pixel has a density value)
      for( G4int ww = 0; ww < rows; ww++) {
	for( G4int xx = 0; xx < columns; xx++) {
	  mean = tab[ww][xx];
	  density = Pixel2density(mean);
	  std::fwrite(&density, sizeof(G4float), 1, fileOut);
	}
      }
      
    } else {
      // density value is the average of a square region of
      // compression*compression pixels
      for(G4int ww = 0; ww < rows ;ww += compSize ) {
	for(G4int xx = 0; xx < columns ;xx +=compSize ) {
	  overflow = false;
	  mean = 0;
	  for(int sumx = 0; sumx < compSize; sumx++) {
	    for(int sumy = 0; sumy < compSize; sumy++) {
	      if(ww+sumy >= rows || xx+sumx >= columns) overflow = true;
	      mean += tab[ww+sumy][xx+sumx];
	    }
	    if(overflow) break;
	  }
	  mean /= compSize*compSize;
	  cpt = 1;
	  
	  if(!overflow) {
	    density = Pixel2density(mean);
	    std::fwrite(&density, sizeof(G4float), 1, fileOut);
	  }
	}
	
      }
    }
    
    std::fclose(fileOut);
    
    delete [] nameProcessed;

    /*    for ( G4int i = 0; i < rows; i ++ ) {
      delete [] tab[i];
    }
    delete [] tab;
    */

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

G4float DicomHandler::Pixel2density(G4int pixel)
{
    G4float density = -1.;
    G4int nbrequali = 0;
    G4double deltaCT = 0;
    G4double deltaDensity = 0;

    // CT2Density.dat contains the calibration curve to convert CT (Hounsfield)
    // number to physical density
    std::ifstream calibration(this->calibrationDensityFileName);
    calibration >> nbrequali;

    G4double * valuedensity = new G4double[nbrequali];
    G4double * valueCT = new G4double[nbrequali];

    if(!calibration) {
	G4cerr << "@@@ No value to transform pixels in density!" << G4endl;
	exit(1);

    } else { // calibration was successfully opened
	for(G4int i = 0; i < nbrequali; i++) { // Loop to store all the pts in CT2Density.dat
	    calibration >> valueCT[i] >> valuedensity[i];
	}
    }
    calibration.close();

    for(G4int j = 1; j < nbrequali; j++) {
	if( pixel >= valueCT[j-1] && pixel < valueCT[j]) {

	    deltaCT = valueCT[j] - valueCT[j-1];
	    deltaDensity = valuedensity[j] - valuedensity[j-1];

	    // interpolating linearly
	    density = valuedensity[j] - ((valueCT[j] - pixel)*deltaDensity/deltaCT );
	    break;
	}
    }

    if(density < 0.) {
	std::printf("@@@ Error density = %f && Pixel = %i (0x%x) && deltaDensity/deltaCT = %f\n",density,pixel,pixel, deltaDensity/deltaCT);
    }
    
    delete [] valuedensity;
    delete [] valueCT;

    return density;
}


void DicomHandler::CheckFileFormat(G4String directory, G4String fileName, G4String calibrationDensityFileName)
{
	this->dicomDirectory=directory;
        this->calibrationDensityFileName=this->dicomDirectory+calibrationDensityFileName;
        G4String fullFileName=this->dicomDirectory+fileName;

	std::cout << "fullFileName "<< fullFileName<< G4endl;

	char * fullName = new char[FILENAMESIZE];// added the full path
	std::ifstream checkData(fullFileName);
    char * oneLine = new char[128];

    if(!(checkData.is_open())) { //Check existance of Data.dat

	G4cout << "\nDicomG4 needs Data.dat :\n\tFirst line: number of image pixel for a "
	       << "voxel (G4Box)\n\tSecond line: number of images (CT slices) to "
	       << "read\n\tEach following line contains the name of a Dicom image except "
	       << "for the .dcm extension\n";
	exit(0);
    }

    checkData >> compression;
    checkData >> nFiles;
    G4String oneName;
    checkData.getline(oneLine,100);
    std::ifstream testExistence;
    G4bool existAlready = true;
    for(G4int rep = 0; rep < nFiles; rep++) { 
      checkData.getline(oneLine,100);
      oneName = oneLine;
      oneName += ".g4dcm"; // create dicomFile.g4dcm
      G4cout << nFiles << " test file " << oneName << G4endl;
//      testExistence.open(oneName.data()); // added the full path
	strcpy(fullName,this->dicomDirectory);	
	strcat(fullName,oneName.data());
      testExistence.open(fullName);  // added the full path
      if(!(testExistence.is_open())) {
	existAlready = false;
	testExistence.clear();
	testExistence.close();
      }
      testExistence.clear();
      testExistence.close();
    }

    ReadMaterialIndices( checkData );

    checkData.close();
    delete [] oneLine;

   if( existAlready == false ) { // The files *.g4dcm have to be created

 	G4cout << "\nAll the necessary images were not found in processed form, starting "
	       << "with .dcm images\n";

	FILE * dicom;
	FILE * lecturePref;
	char * compressionc = new char[LINEBUFFSIZE];
	char * maxc = new char[LINEBUFFSIZE];
	//char name[300], inputFile[300];
	char * name = new char[FILENAMESIZE];
	char * inputFile = new char[FILENAMESIZE];


// 	lecturePref = std::fopen("Data.dat","r"); // added the full path
	lecturePref = std::fopen(fullFileName,"r"); // added the full path
	std::fscanf(lecturePref,"%s",compressionc);
	compression = atoi(compressionc);
	std::fscanf(lecturePref,"%s",maxc);
	nFiles = atoi(maxc);
	G4cout << " nFiles " << nFiles << G4endl;
	
	for( G4int i = 1; i <= nFiles; i++ ) { // Begin loop on filenames

	    std::fscanf(lecturePref,"%s",inputFile);
	    std::sprintf(name,"%s.dcm",inputFile);
	    std::cout << "check 1: " << name << std::endl;
	    //  Open input file and give it to gestion_dicom :
	    std::printf("### Opening %s and reading :\n",name);
// 	    dicom = std::fopen(name,"rb"); // added the full path
	    strcpy(fullName,this->dicomDirectory);	
	    strcat(fullName,name);
	    dicom = std::fopen(fullName,"rb"); // added the full path
	    // Reading the .dcm in two steps:
	    //      1.  reading the header
	    //	2. reading the pixel data and store the density in Moyenne.dat
	    if( dicom != 0 ) {
		ReadFile(dicom,inputFile);
	    } else {
		G4cout << "\nError opening file : " << name << G4endl;
	    }
	    std::fclose(dicom);
	}
	std::fclose(lecturePref);

	delete [] compressionc;
	delete [] maxc;
	delete [] name;
	delete [] inputFile;

    } 

}


template <class Type>
void DicomHandler::GetValue(char * _val, Type & _rval) {

#if BYTE_ORDER == BIG_ENDIAN
    if(littleEndian) {      // little endian
#else // BYTE_ORDER == LITTLE_ENDIAN
    if(!littleEndian) {     // big endian
#endif
	const int SIZE = sizeof(_rval);
	char ctemp;
	for(int i = 0; i < SIZE/2; i++) {
	    ctemp = _val[i];
	    _val[i] = _val[SIZE - 1 - i];
	    _val[SIZE - 1 - i] = ctemp;
	}
    }
    _rval = *(Type *)_val;
}
