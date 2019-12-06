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
//
/// \file medical/DICOM/src/DicomHandler.cc
/// \brief Implementation of the DicomHandler class
//
// The code was written by :
//      *Louis Archambault louis.archambault@phy.ulaval.ca,
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
// + University Laval, Quebec (QC) Canada
//*******************************************************
//
//*******************************************************
//
/// DicomHandler.cc :
///        - Handling of DICM images
///         - Reading headers and pixels
///        - Transforming pixel to density and creating *.g4dcm
///          files
//*******************************************************

#include "DicomHandler.hh"
#include "globals.hh"
#include "G4ios.hh"
#include <fstream>

#include <cctype>
#include <cstring>

#include "DicomPhantomZSliceHeader.hh"
#include "DicomPhantomZSliceMerged.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DicomHandler* DicomHandler::fInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DicomHandler* DicomHandler::Instance()
{
  if (fInstance == 0)
  {
    static DicomHandler dicomhandler;
    fInstance = &dicomhandler;
  }
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String DicomHandler::GetDicomDataPath()
{
    // default is current directory
    G4String driverPath = ".";
    // check environment
    const char* path = std::getenv("DICOM_PATH");

    if(path)
    {
        // if is set in environment
        return G4String(path);
    }
    else
    {
        // if DICOM_USE_HEAD, look for data installation
#ifdef DICOM_USE_HEAD
        G4cerr << "Warning! DICOM was compiled with DICOM_USE_HEAD option but "
               << "the DICOM_PATH was not set!" << G4endl;
        G4String datadir = G4GetEnv<G4String>("G4ENSDFSTATEDATA", "");
        if(datadir.length() > 0)
        {
            auto _last = datadir.rfind("/");
            if(_last != std::string::npos)
                datadir.erase(_last);
            driverPath = datadir + "/DICOM1.1/DICOM_HEAD_TEST";
            G4int rc = setenv("DICOM_PATH", driverPath.c_str(), 0);
            G4cerr << "\t --> Using '" << driverPath << "'..." << G4endl;
            G4ConsumeParameters(rc);
        }
#endif
    }
    return driverPath;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String DicomHandler::GetDicomDataFile()
{
#if defined(DICOM_USE_HEAD) && defined(G4_DCMTK)
    return GetDicomDataPath() + "/Data.dat.new";
#else
    return GetDicomDataPath() + "/Data.dat";
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef DICOM_USE_HEAD
DicomHandler::DicomHandler()
:DATABUFFSIZE(8192), LINEBUFFSIZE(5020), FILENAMESIZE(512),
 fCompression(0),fNFiles(0), fRows(0), fColumns(0), 
 fBitAllocated(0), fMaxPixelValue(0),
 fMinPixelValue(0),fPixelSpacingX(0.),
 fPixelSpacingY(0.),fSliceThickness(0.), 
 fSliceLocation(0.),fRescaleIntercept(0), 
 fRescaleSlope(0),fLittleEndian(true), 
 fImplicitEndian(false),fPixelRepresentation(0), 
 fNbrequali(0),fValueDensity(NULL),fValueCT(NULL),fReadCalibration(false),
 fMergedSlices(NULL),fCt2DensityFile("null.dat")
{
    fMergedSlices = new DicomPhantomZSliceMerged;
    fDriverFile = GetDicomDataFile();
    G4cout << "Reading the DICOM_HEAD project " << fDriverFile << G4endl;
}
#else
DicomHandler::DicomHandler()
:DATABUFFSIZE(8192), LINEBUFFSIZE(5020), FILENAMESIZE(512),
 fCompression(0), fNFiles(0), fRows(0), fColumns(0),
 fBitAllocated(0), fMaxPixelValue(0), fMinPixelValue(0),
 fPixelSpacingX(0.), fPixelSpacingY(0.),fSliceThickness(0.), 
 fSliceLocation(0.), fRescaleIntercept(0), fRescaleSlope(0),
 fLittleEndian(true),fImplicitEndian(false),fPixelRepresentation(0),
 fNbrequali(0),fValueDensity(NULL),fValueCT(NULL),
 fReadCalibration(false),fMergedSlices(NULL),
 fDriverFile("Data.dat"),fCt2DensityFile("CT2Density.dat")
{
 fMergedSlices = new DicomPhantomZSliceMerged;
}
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DicomHandler::~DicomHandler()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DicomHandler::ReadFile(FILE* dicom, char* filename2)
{
    G4cout << " ReadFile " << filename2 << G4endl;

    G4int returnvalue = 0; size_t rflag = 0;
    char * buffer = new char[LINEBUFFSIZE];

    fImplicitEndian = false;
    fLittleEndian = true;

    rflag = std::fread( buffer, 1, 128, dicom ); // The first 128 bytes
                                                 //are not important
                                                 // Reads the "DICOM" letters
    rflag = std::fread( buffer, 1, 4, dicom );
    // if there is no preamble, the FILE pointer is rewinded.
    if(std::strncmp("DICM", buffer, 4) != 0) {
        std::fseek(dicom, 0, SEEK_SET);
        fImplicitEndian = true;
    }

    short readGroupId;    // identify the kind of input data
    short readElementId;  // identify a particular type information
    short elementLength2; // deal with element length in 2 bytes
                          //unsigned int elementLength4; 
    // deal with element length in 4 bytes
    unsigned long elementLength4; // deal with element length in 4 bytes

    char * data = new char[DATABUFFSIZE];

    // Read information up to the pixel data
    while(true) {
      
      //Reading groups and elements :
      readGroupId = 0;
      readElementId = 0;
      // group ID
      rflag = std::fread(buffer, 2, 1, dicom);
      GetValue(buffer, readGroupId);
      // element ID
      rflag = std::fread(buffer, 2, 1, dicom);
      GetValue(buffer, readElementId);
      
      // Creating a tag to be identified afterward
      G4int tagDictionary = readGroupId*0x10000 + readElementId;
      
      // beginning of the pixels
      if(tagDictionary == 0x7FE00010) {
        // Following 2 fread's are modifications to original DICOM example
        // (Jonathan Madsen)
        if(!fImplicitEndian)
            rflag = std::fread(buffer,2,1,dicom);   // Reserved 2 bytes
        // (not used for pixels)
        rflag = std::fread(buffer,4,1,dicom);   // Element Length  
        // (not used for pixels)
        break;      // Exit to ReadImageData()
      }
      
      // VR or element length
      rflag = std::fread(buffer,2,1,dicom);
      GetValue(buffer, elementLength2);
      
      // If value representation (VR) is OB, OW, SQ, UN, added OF and UT
      //the next length is 32 bits
      if((elementLength2 == 0x424f ||     // "OB"
        elementLength2 == 0x574f ||     // "OW"
        elementLength2 == 0x464f ||     // "OF"
        elementLength2 == 0x5455 ||     // "UT"
        elementLength2 == 0x5153 ||     // "SQ"
        elementLength2 == 0x4e55) &&    // "UN"
       !fImplicitEndian ) {             // explicit VR
      
      rflag = std::fread(buffer, 2, 1, dicom); // Skip 2 reserved bytes
      
      // element length
      rflag = std::fread(buffer, 4, 1, dicom);
      GetValue(buffer, elementLength4);
      
      if(elementLength2 == 0x5153)
        {
          if(elementLength4 == 0xFFFFFFFF)
            {
            read_undefined_nested( dicom );
            elementLength4=0;
            }  else{
            if(read_defined_nested( dicom, elementLength4 )==0){
            G4Exception("DicomHandler::ReadFile",
                      "DICOM001",
                      FatalException,
                      "Function read_defined_nested() failed!");
            }
          }
        } else  {
        // Reading the information with data
        rflag = std::fread(data, elementLength4,1,dicom);
      }
      
      }  else {

      //  explicit with VR different than previous ones
      if(!fImplicitEndian || readGroupId == 2) {
        
        //G4cout << "Reading  DICOM files with Explicit VR"<< G4endl;
        // element length (2 bytes)
        rflag = std::fread(buffer, 2, 1, dicom);
        GetValue(buffer, elementLength2);
        elementLength4 = elementLength2;
        
        rflag = std::fread(data, elementLength4, 1, dicom);
        
      } else {                                  // Implicit VR
        
        //G4cout << "Reading  DICOM files with Implicit VR"<< G4endl;
        
        // element length (4 bytes)
        if(std::fseek(dicom, -2, SEEK_CUR) != 0) {
          G4Exception("DicomHandler::ReadFile",
                  "DICOM001",
                  FatalException,
                  "fseek failed");
        }
    
        rflag = std::fread(buffer, 4, 1, dicom);
        GetValue(buffer, elementLength4);
        
        //G4cout <<  std::hex<< elementLength4 << G4endl;
        
        if(elementLength4 == 0xFFFFFFFF)
          {
            read_undefined_nested(dicom);
            elementLength4=0;
          }  else{
          rflag = std::fread(data, elementLength4, 1, dicom);
        }
        
      }
      }
      
      // NULL termination
      data[elementLength4] = '\0';
      
      // analyzing information
      GetInformation(tagDictionary, data);
    }
    
    G4String fnameG4DCM = G4String(filename2) + ".g4dcm";

    // Perform functions originally written straight to file
    DicomPhantomZSliceHeader* zslice = new DicomPhantomZSliceHeader(fnameG4DCM);
    for(auto ite = fMaterialIndices.cbegin();
             ite != fMaterialIndices.cend(); ++ite)
    {
      zslice->AddMaterial(ite->second);
    }
    
    zslice->SetNoVoxelX(fColumns/fCompression);
    zslice->SetNoVoxelY(fRows/fCompression);
    zslice->SetNoVoxelZ(1);

    zslice->SetMinX(-fPixelSpacingX*fColumns/2.);
    zslice->SetMaxX(fPixelSpacingX*fColumns/2.);

    zslice->SetMinY(-fPixelSpacingY*fRows/2.);
    zslice->SetMaxY(fPixelSpacingY*fRows/2.);

    zslice->SetMinZ(fSliceLocation-fSliceThickness/2.);
    zslice->SetMaxZ(fSliceLocation+fSliceThickness/2.);

    //===

    ReadData( dicom, filename2 );

    // DEPRECIATED
    //StoreData( foutG4DCM );
    //foutG4DCM.close();

    StoreData( zslice );

    // Dumped 2 file after DicomPhantomZSliceMerged has checked for consistency
    //zslice->DumpToFile();

    fMergedSlices->AddZSlice(zslice);

    //
    delete [] buffer;
    delete [] data;

    if (rflag) return returnvalue;
    return returnvalue;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DicomHandler::GetInformation(G4int & tagDictionary, char * data)
{
    if(tagDictionary == 0x00280010 ) { // Number of Rows
        GetValue(data, fRows);
        std::printf("[0x00280010] Rows -> %i\n",fRows);

    } else if(tagDictionary == 0x00280011 ) { // Number of fColumns
        GetValue(data, fColumns);
        std::printf("[0x00280011] Columns -> %i\n",fColumns);

    } else if(tagDictionary == 0x00280102 ) { // High bits  ( not used )
        short highBits;
        GetValue(data, highBits);
        std::printf("[0x00280102] High bits -> %i\n",highBits);

    } else if(tagDictionary == 0x00280100 ) { // Bits allocated
        GetValue(data, fBitAllocated);
        std::printf("[0x00280100] Bits allocated -> %i\n", fBitAllocated);

    } else if(tagDictionary == 0x00280101 ) { //  Bits stored ( not used )
        short bitStored;
        GetValue(data, bitStored);
        std::printf("[0x00280101] Bits stored -> %i\n",bitStored);

    } else if(tagDictionary == 0x00280106 ) { //  Min. pixel value
        GetValue(data, fMinPixelValue);
        std::printf("[0x00280106] Min. pixel value -> %i\n", fMinPixelValue);

    } else if(tagDictionary == 0x00280107 ) { //  Max. pixel value
        GetValue(data, fMaxPixelValue);
        std::printf("[0x00280107] Max. pixel value -> %i\n", fMaxPixelValue);

    } else if(tagDictionary == 0x00281053) { //  Rescale slope
        fRescaleSlope = atoi(data);
        std::printf("[0x00281053] Rescale Slope -> %d\n", fRescaleSlope);

    } else if(tagDictionary == 0x00281052 ) { // Rescalse intercept
        fRescaleIntercept = atoi(data);
        std::printf("[0x00281052] Rescale Intercept -> %d\n", 
                    fRescaleIntercept );

    } else if(tagDictionary == 0x00280103 ) {
        //  Pixel representation ( functions not design to read signed bits )
      fPixelRepresentation = atoi(data); // 0: unsigned  1: signed
        std::printf("[0x00280103] Pixel Representation -> %i\n",
                    fPixelRepresentation);
        if(fPixelRepresentation == 1 ) {
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

    } else if(tagDictionary == 0x00283000 ) { //Modality LUT Sequence(not used)
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

    } else if(tagDictionary == 0x00280120 ) { // Pixel Padding Value (not used)
      std::printf("[0x00280120] Pixel Padding Value US or SS 1 -> %s\n", data);

    } else if(tagDictionary == 0x00280030 ) { // Pixel Spacing
        G4String datas(data);
        G4int iss = G4int(datas.find('\\'));
        fPixelSpacingX = atof( datas.substr(0,iss).c_str() );
        fPixelSpacingY = atof( datas.substr(iss+1,datas.length()).c_str() );

    } else if(tagDictionary == 0x00200037 ) { // Image Orientation ( not used )
        std::printf("[0x00200037] Image Orientation (Phantom) -> %s\n", data);

    } else if(tagDictionary == 0x00200032 ) { // Image Position ( not used )
        std::printf("[0x00200032] Image Position (Phantom,mm) -> %s\n", data);

    } else if(tagDictionary == 0x00180050 ) { // Slice Thickness
      fSliceThickness = atof(data);
      std::printf("[0x00180050] Slice Thickness (mm) -> %f\n", fSliceThickness);

    } else if(tagDictionary == 0x00201041 ) { // Slice Location
      fSliceLocation = atof(data);
      std::printf("[0x00201041] Slice Location -> %f\n", fSliceLocation);

    } else if(tagDictionary == 0x00280004 ) { // Photometric Interpretation
      // ( not used )
        std::printf("[0x00280004] Photometric Interpretation -> %s\n", data);

    } else if(tagDictionary == 0x00020010) { // Endian
        if(strcmp(data, "1.2.840.10008.1.2") == 0)
            fImplicitEndian = true;
        else if(strncmp(data, "1.2.840.10008.1.2.2", 19) == 0)
            fLittleEndian = false;
        //else 1.2.840..10008.1.2.1 (explicit little endian)

        std::printf("[0x00020010] Endian -> %s\n", data);
    }

    // others
    else {
        //std::printf("[0x%x] -> %s\n", tagDictionary, data);
        ;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DicomHandler::StoreData(DicomPhantomZSliceHeader* dcmPZSH)
{
    G4int mean;
    G4double density;
    G4bool overflow = false;

    if(!dcmPZSH) { return; }

    dcmPZSH->SetSliceLocation(fSliceLocation);

    //----- Print indices of material
    if(fCompression == 1) { // no fCompression: each pixel has a density value)
        for( G4int ww = 0; ww < fRows; ++ww) {
            dcmPZSH->AddRow();
            for( G4int xx = 0; xx < fColumns; ++xx) {
                mean = fTab[ww][xx];
                density = Pixel2density(mean);
                dcmPZSH->AddValue(density);
                dcmPZSH->AddMateID(GetMaterialIndex(G4float(density)));
            }
        }

    } else {
        // density value is the average of a square region of
        // fCompression*fCompression pixels
      for(G4int ww = 0; ww < fRows ;ww += fCompression ) {
        dcmPZSH->AddRow();
        for(G4int xx = 0; xx < fColumns ;xx +=fCompression ) {
          overflow = false;
          mean = 0;
          for(G4int sumx = 0; sumx < fCompression; ++sumx) {
            for(G4int sumy = 0; sumy < fCompression; ++sumy) {
              if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
              mean += fTab[ww+sumy][xx+sumx];
            }
            if(overflow) break;
          }
          mean /= fCompression*fCompression;
          
          if(!overflow) {
            density = Pixel2density(mean);
            dcmPZSH->AddValue(density);
            dcmPZSH->AddMateID(GetMaterialIndex(G4float(density)));
          }
        }
      }
    }
    
    dcmPZSH->FlipData();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// This function is depreciated as it is handled by 
// DicomPhantomZSliceHeader::DumpToFile
void DicomHandler::StoreData(std::ofstream& foutG4DCM)
{
  G4int mean;
  G4double density;
  G4bool overflow = false;
  
  //----- Print indices of material
  if(fCompression == 1) { // no fCompression: each pixel has a density value)
    for( G4int ww = 0; ww < fRows; ++ww) {
      for( G4int xx = 0; xx < fColumns; ++xx) {
        mean = fTab[ww][xx];
        density = Pixel2density(mean);
        foutG4DCM << GetMaterialIndex( G4float(density) ) << " ";
      }
      foutG4DCM << G4endl;
    }
    
  } else {
    // density value is the average of a square region of
    // fCompression*fCompression pixels
    for(G4int ww = 0; ww < fRows ;ww += fCompression ) {
      for(G4int xx = 0; xx < fColumns ;xx +=fCompression ) {
        overflow = false;
        mean = 0;
        for(G4int sumx = 0; sumx < fCompression; ++sumx) {
          for(G4int sumy = 0; sumy < fCompression; ++sumy) {
            if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
            mean += fTab[ww+sumy][xx+sumx];
          }
          if(overflow) break;
        }
        mean /= fCompression*fCompression;
        
        if(!overflow) {
          density = Pixel2density(mean);
          foutG4DCM << GetMaterialIndex( G4float(density) ) << " ";
        }
      }
      foutG4DCM << G4endl;
    }
    
  }
  
  //----- Print densities
  if(fCompression == 1) { // no fCompression: each pixel has a density value)
    for( G4int ww = 0; ww < fRows; ww++) {
      for( G4int xx = 0; xx < fColumns; xx++) {
        mean = fTab[ww][xx];
        density = Pixel2density(mean);
        foutG4DCM << density << " ";
        if( xx%8 == 3 ) foutG4DCM << G4endl; // just for nicer reading
      }
    }
    
  } else {
    // density value is the average of a square region of
    // fCompression*fCompression pixels
    for(G4int ww = 0; ww < fRows ;ww += fCompression ) {
      for(G4int xx = 0; xx < fColumns ;xx +=fCompression ) {
        overflow = false;
        mean = 0;
        for(G4int sumx = 0; sumx < fCompression; ++sumx) {
          for(G4int sumy = 0; sumy < fCompression; ++sumy) {
            if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
            mean += fTab[ww+sumy][xx+sumx];
          }
          if(overflow) break;
        }
        mean /= fCompression*fCompression;
        
        if(!overflow) {
          density = Pixel2density(mean);
          foutG4DCM << density  << " ";
          if( xx/fCompression%8 == 3 ) foutG4DCM << G4endl; // just for nicer
          // reading
        }
      }
    }
    
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
void DicomHandler::ReadMaterialIndices( std::ifstream& finData)
{
  unsigned int nMate;
  G4String mateName;
  G4float densityMax;
  finData >> nMate;
  if( finData.eof() ) return;
  
  G4cout << " ReadMaterialIndices " << nMate << G4endl;
  for( unsigned int ii = 0; ii < nMate; ++ii )
  {
    finData >> mateName >> densityMax;
    fMaterialIndices[densityMax] = mateName;
    //    G4cout << ii << " ReadMaterialIndices " << mateName << " " 
    //     << densityMax << G4endl;
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unsigned int DicomHandler::GetMaterialIndex( G4float density )
{
  std::map<G4float,G4String>::const_reverse_iterator ite;
  std::size_t ii = fMaterialIndices.size();
 
  for( ite = fMaterialIndices.crbegin(); ite != fMaterialIndices.crend();
       ++ite, ii-- )
  {
    if( density >= (*ite).first )  { break; }
  }
  
  if(ii == fMaterialIndices.size())  { ii = fMaterialIndices.size()-1; }
  
  return unsigned(ii);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DicomHandler::ReadData(FILE *dicom,char * filename2)
{
  G4int returnvalue = 0; size_t rflag = 0;
  
  //  READING THE PIXELS :
  G4int w = 0;
  
  fTab = new G4int*[fRows];
  for ( G4int i = 0; i < fRows; ++i )
  {
    fTab[i] = new G4int[fColumns];
  }
  
  if(fBitAllocated == 8) { // Case 8 bits :
    
    std::printf("@@@ Error! Picture != 16 bits...\n");
    std::printf("@@@ Error! Picture != 16 bits...\n");
    std::printf("@@@ Error! Picture != 16 bits...\n");
    
    unsigned char ch = 0;
    
    for(G4int j = 0; j < fRows; ++j) {
      for(G4int i = 0; i < fColumns; ++i) {
        w++;
        rflag = std::fread( &ch, 1, 1, dicom);
        fTab[j][i] = ch*fRescaleSlope + fRescaleIntercept;
      }
    }
    returnvalue = 1;
    
  } else { //  from 12 to 16 bits :
    char sbuff[2];
    short pixel;
    for( G4int j = 0; j < fRows; ++j) {
      for( G4int i = 0; i < fColumns; ++i) {
        w++;
        rflag = std::fread(sbuff, 2, 1, dicom);
        GetValue(sbuff, pixel);
        fTab[j][i] = pixel*fRescaleSlope + fRescaleIntercept;
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
  rflag = std::fwrite(&nMate, sizeof(unsigned int), 1, fileOut);
  //--- Write materials
  for( auto ite = fMaterialIndices.cbegin();
            ite != fMaterialIndices.cend(); ++ite )
  {
    G4String mateName = (*ite).second;
    for( std::size_t ii = (*ite).second.length(); ii < 40; ++ii )
    {
      mateName += " ";
    }         //mateName = const_cast<char*>(((*ite).second).c_str());
    
    const char* mateNameC = mateName.c_str();
    rflag = std::fwrite(mateNameC, sizeof(char),40, fileOut);
  }

  unsigned int fRowsC = fRows/fCompression;
  unsigned int fColumnsC = fColumns/fCompression;
  unsigned int planesC = 1;
  G4float pixelLocationXM = -G4float(fPixelSpacingX*fColumns/2.);
  G4float pixelLocationXP = G4float(fPixelSpacingX*fColumns/2.);
  G4float pixelLocationYM = -G4float(fPixelSpacingY*fRows/2.);
  G4float pixelLocationYP = G4float(fPixelSpacingY*fRows/2.);
  G4float fSliceLocationZM = G4float(fSliceLocation-fSliceThickness/2.);
  G4float fSliceLocationZP = G4float(fSliceLocation+fSliceThickness/2.);
  //--- Write number of voxels (assume only one voxel in Z)
  rflag = std::fwrite(&fRowsC, sizeof(unsigned int), 1, fileOut);
  rflag = std::fwrite(&fColumnsC, sizeof(unsigned int), 1, fileOut);
  rflag = std::fwrite(&planesC, sizeof(unsigned int), 1, fileOut);
  //--- Write minimum and maximum extensions
  rflag = std::fwrite(&pixelLocationXM, sizeof(G4float), 1, fileOut);
  rflag = std::fwrite(&pixelLocationXP, sizeof(G4float), 1, fileOut);
  rflag = std::fwrite(&pixelLocationYM, sizeof(G4float), 1, fileOut);
  rflag = std::fwrite(&pixelLocationYP, sizeof(G4float), 1, fileOut);
  rflag = std::fwrite(&fSliceLocationZM, sizeof(G4float), 1, fileOut);
  rflag = std::fwrite(&fSliceLocationZP, sizeof(G4float), 1, fileOut);
  // rflag = std::fwrite(&fCompression, sizeof(unsigned int), 1, fileOut);
  
  std::printf("%8i   %8i\n",fRows,fColumns);
  std::printf("%8f   %8f\n",fPixelSpacingX,fPixelSpacingY);
  std::printf("%8f\n", fSliceThickness);
  std::printf("%8f\n", fSliceLocation);
  std::printf("%8i\n", fCompression);
  
  G4int compSize = fCompression;
  G4int mean;
  G4float density;
  G4bool overflow = false;
  
  //----- Write index of material for each pixel
  if(compSize == 1) { // no fCompression: each pixel has a density value)
    for( G4int ww = 0; ww < fRows; ++ww) {
      for( G4int xx = 0; xx < fColumns; ++xx) {
        mean = fTab[ww][xx];
        density = Pixel2density(mean);
        unsigned int mateID = GetMaterialIndex( density );
        rflag = std::fwrite(&mateID, sizeof(unsigned int), 1, fileOut);
      }
    }
    
  } else {
    // density value is the average of a square region of
    // fCompression*fCompression pixels
    for(G4int ww = 0; ww < fRows ;ww += compSize ) {
      for(G4int xx = 0; xx < fColumns ;xx +=compSize ) {
        overflow = false;
        mean = 0;
        for(G4int sumx = 0; sumx < compSize; ++sumx) {
          for(G4int sumy = 0; sumy < compSize; ++sumy) {
            if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
            mean += fTab[ww+sumy][xx+sumx];
          }
          if(overflow) break;
        }
        mean /= compSize*compSize;
        
        if(!overflow) {
          density = Pixel2density(mean);
          unsigned int mateID = GetMaterialIndex( density );
          rflag = std::fwrite(&mateID, sizeof(unsigned int), 1, fileOut);
        }
      }
      
    }
  }
  
  //----- Write density for each pixel
  if(compSize == 1) { // no fCompression: each pixel has a density value)
    for( G4int ww = 0; ww < fRows; ++ww) {
      for( G4int xx = 0; xx < fColumns; ++xx) {
        mean = fTab[ww][xx];
        density = Pixel2density(mean);
        rflag = std::fwrite(&density, sizeof(G4float), 1, fileOut);
      }
    }
    
  } else {
    // density value is the average of a square region of
    // fCompression*fCompression pixels
    for(G4int ww = 0; ww < fRows ;ww += compSize ) {
      for(G4int xx = 0; xx < fColumns ;xx +=compSize ) {
        overflow = false;
        mean = 0;
        for(G4int sumx = 0; sumx < compSize; ++sumx) {
          for(G4int sumy = 0; sumy < compSize; ++sumy) {
            if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
            mean += fTab[ww+sumy][xx+sumx];
          }
          if(overflow) break;
        }
        mean /= compSize*compSize;
        
        if(!overflow) {
          density = Pixel2density(mean);
          rflag = std::fwrite(&density, sizeof(G4float), 1, fileOut);
        }
      }
      
    }
  }
  
  rflag = std::fclose(fileOut);
  
  delete [] nameProcessed;

  /*    for ( G4int i = 0; i < fRows; i ++ ) {
        delete [] fTab[i];
        }
     delete [] fTab;
  */
  
  if (rflag) return returnvalue;
  return returnvalue;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.......

// DICOM HEAD does not use the calibration curve

#ifdef DICOM_USE_HEAD
void DicomHandler::ReadCalibration()
{
fNbrequali = 0;
fReadCalibration = false;
G4cout << "No calibration curve for the DICOM HEAD needed!" << G4endl;
}
#else
// Separated out of Pixel2density
// No need to read in same calibration EVERY time
// Increases the speed of reading file by several orders of magnitude

void DicomHandler::ReadCalibration()
{
 fNbrequali = 0;
// CT2Density.dat contains the calibration curve to convert CT (Hounsfield)
// number to physical density
 std::ifstream calibration(fCt2DensityFile.c_str());
 calibration >> fNbrequali; 
 fValueDensity = new G4double[fNbrequali];
 fValueCT = new G4double[fNbrequali];

 if(!calibration) {
   G4Exception("DicomHandler::ReadFile","DICOM001", FatalException,
               "@@@ No value to transform pixels in density!");
   } 
   else { // calibration was successfully opened
      for(G4int i = 0; i < fNbrequali; ++i) 
         { // Loop to store all the pts in CT2Density.dat
        calibration >> fValueCT[i] >> fValueDensity[i];
        }
  }
calibration.close();
fReadCalibration = true;
}
#endif

#ifdef DICOM_USE_HEAD
G4float DicomHandler::Pixel2density(G4int pixel)
{
 G4float density = -1;

//Air
  if (pixel == -998.) density = 0.001290;
//Soft Tissue
  else if ( pixel == 24.) density = 1.00;
//Brain
  else if ( pixel == 52.) density = 1.03; 
// Spinal disc
   else if ( pixel == 92.) density = 1.10;
// Trabecular bone
   else if ( pixel == 197.) density = 1.18;
// Cortical Bone
   else if ( pixel == 923.) density = 1.85;
// Tooth dentine
   else if ( pixel == 1280.) density = 2.14;
//Tooth enamel
   else if ( pixel == 2310.) density = 2.89; 

return density;
}

#else
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4float DicomHandler::Pixel2density(G4int pixel)
{
  if(!fReadCalibration) { ReadCalibration(); }
  
  G4float density = -1.;
  G4double deltaCT = 0;
  G4double deltaDensity = 0;
  
  
  for(G4int j = 1; j < fNbrequali; ++j) {
    if( pixel >= fValueCT[j-1] && pixel < fValueCT[j]) {
      
      deltaCT = fValueCT[j] - fValueCT[j-1];
      deltaDensity = fValueDensity[j] - fValueDensity[j-1];
      
      // interpolating linearly
      density = G4float(fValueDensity[j]
                        -((fValueCT[j] - pixel)*deltaDensity/deltaCT ));
      break;
    }
  }
  
  if(density < 0.) {
    std::printf("@@@ Error density = %f && Pixel = %i \
      (0x%x) && deltaDensity/deltaCT = %f\n",density,pixel,pixel,
                deltaDensity/deltaCT);
  }
  
  return density;
}
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DicomHandler::CheckFileFormat()
{
  std::ifstream checkData(fDriverFile.c_str());
  char * oneLine = new char[128];
  
  if(!(checkData.is_open())) { //Check existance of Data.dat
    
    G4String message = 
      "\nDicomG4 needs Data.dat (or another driver file specified";
    message += " in command line):\n";
    message += "\tFirst line: number of image pixel for a voxel (G4Box)\n";
    message += "\tSecond line: number of images (CT slices) to read\n";
    message += "\tEach following line contains the name of a Dicom image";
    message += " except for the .dcm extension";
    G4Exception("DicomHandler::ReadFile",
                "DICOM001",
                FatalException,
                message.c_str());
  }
  
  checkData >> fCompression;
  checkData >> fNFiles;
  G4String oneName;
  checkData.getline(oneLine,100);
  std::ifstream testExistence;
  G4bool existAlready = true;
  for(G4int rep = 0; rep < fNFiles; ++rep) {
    checkData.getline(oneLine,100);
    oneName = oneLine;
    oneName += ".g4dcm"; // create dicomFile.g4dcm
    G4cout << fNFiles << " test file " << oneName << G4endl;
    testExistence.open(oneName.data());
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
  
  if( existAlready == false  ) { // The files *.g4dcm have to be created
    
    G4cout << "\nAll the necessary images were not found in processed form "
           << ", starting with .dcm images\n";
    
    FILE * dicom;
    FILE * lecturePref;
    char * fCompressionc = new char[LINEBUFFSIZE];
    char * maxc = new char[LINEBUFFSIZE];
    //char name[300], inputFile[300];
    char * name = new char[FILENAMESIZE];
    char * inputFile = new char[FILENAMESIZE];
    G4int rflag;
    lecturePref = std::fopen(fDriverFile.c_str(),"r");
 
    rflag = std::fscanf(lecturePref,"%s",fCompressionc);
    fCompression = atoi(fCompressionc);
    rflag = std::fscanf(lecturePref,"%s",maxc);
    fNFiles = atoi(maxc);
    G4cout << " fNFiles " << fNFiles << G4endl;
  
///////////////////// 

#ifdef DICOM_USE_HEAD
    for( G4int i = 1; i <= fNFiles; ++i ) 
    { // Begin loop on filenames
     rflag = std::fscanf(lecturePref,"%s",inputFile);
      G4String path = GetDicomDataPath() + "/";
    
     std::sprintf(name,"%s.dcm",inputFile);
     //Writes the results to a character string buffer.
     
     char name2[200];
     strcpy(name2, path.c_str());
     strcat(name2, name);
     //  Open input file and give it to gestion_dicom :
     std::printf("### Opening %s and reading :\n",name2);
     dicom = std::fopen(name2,"rb");
     // Reading the .dcm in two steps:
     //      1.  reading the header
     //      2. reading the pixel data and store the density in Moyenne.dat
     if( dicom != 0 ) {
        ReadFile(dicom,inputFile);
      } else {
        G4cout << "\nError opening file : " << name2 << G4endl;
      }
      rflag = std::fclose(dicom);
    }
#else

     for( G4int i = 1; i <= fNFiles; ++i ) 
     { // Begin loop on filenames
      rflag = std::fscanf(lecturePref,"%s",inputFile);
 
      std::sprintf(name,"%s.dcm",inputFile);
      //Writes the results to a character string buffer.
      
       //G4cout << "check: " << name << G4endl;
      //  Open input file and give it to gestion_dicom :
      std::printf("### Opening %s and reading :\n",name);
      dicom = std::fopen(name,"rb");
      // Reading the .dcm in two steps:
      //      1.  reading the header
      //      2. reading the pixel data and store the density in Moyenne.dat
      if( dicom != 0 ) {
        ReadFile(dicom,inputFile);
      } else {
        G4cout << "\nError opening file : " << name << G4endl;
      }
      rflag = std::fclose(dicom);
    }
#endif

    rflag = std::fclose(lecturePref);
    
    // Checks the spacing is correct. Dumps to files
    fMergedSlices->CheckSlices();
    
    delete [] fCompressionc;
    delete [] maxc;
    delete [] name;
    delete [] inputFile;
    if (rflag) return;
    
  }
  
  if(fValueDensity) { delete [] fValueDensity; }
  if(fValueCT) { delete [] fValueCT; }
  if(fMergedSlices) { delete fMergedSlices; }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int DicomHandler::read_defined_nested(FILE * nested,G4int SQ_Length)
{
  //      VARIABLES
  unsigned short item_GroupNumber;
  unsigned short item_ElementNumber;
  G4int item_Length;
  G4int items_array_length=0;
  char * buffer= new char[LINEBUFFSIZE];
  size_t rflag = 0;
  
  while(items_array_length < SQ_Length)
    {
      rflag = std::fread(buffer, 2, 1, nested);
      GetValue(buffer, item_GroupNumber);
      
      rflag = std::fread(buffer, 2, 1, nested);
      GetValue(buffer, item_ElementNumber);
      
      rflag = std::fread(buffer, 4, 1, nested);
      GetValue(buffer, item_Length);
      
      rflag = std::fread(buffer, item_Length, 1, nested);
      
      items_array_length= items_array_length+8+item_Length;
    }
  
  delete [] buffer;
  
  if( SQ_Length>items_array_length )
    return 0;
  else
    return 1;
  if (rflag) return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DicomHandler::read_undefined_nested(FILE * nested)
{
  //      VARIABLES
  unsigned short item_GroupNumber;
  unsigned short item_ElementNumber;
  unsigned int item_Length;
  char * buffer= new char[LINEBUFFSIZE];
  size_t rflag = 0;
  
    do
      {
        rflag = std::fread(buffer, 2, 1, nested);
        GetValue(buffer, item_GroupNumber);
        
        rflag = std::fread(buffer, 2, 1, nested);
        GetValue(buffer, item_ElementNumber);
        
        rflag = std::fread(buffer, 4, 1, nested);
        GetValue(buffer, item_Length);
        
        if(item_Length!=0xffffffff)
          rflag = std::fread(buffer, item_Length, 1, nested);
        else
          read_undefined_item(nested);
        
        
      } while(item_GroupNumber!=0xFFFE || item_ElementNumber!=0xE0DD 
              || item_Length!=0);
    
    delete [] buffer;
    if (rflag) return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DicomHandler::read_undefined_item(FILE * nested)
{
  //      VARIABLES
  unsigned short item_GroupNumber;
  unsigned short item_ElementNumber;
  G4int item_Length; size_t rflag = 0;
  char *buffer= new char[LINEBUFFSIZE];
  
  do
    {
      rflag = std::fread(buffer, 2, 1, nested);
      GetValue(buffer, item_GroupNumber);
      
      rflag = std::fread(buffer, 2, 1, nested);
      GetValue(buffer, item_ElementNumber);
      
      rflag = std::fread(buffer, 4, 1, nested);
      GetValue(buffer, item_Length);
      
      
      if(item_Length!=0)
        rflag = std::fread(buffer,item_Length,1,nested);
      
    }
  while(item_GroupNumber!=0xFFFE || item_ElementNumber!=0xE00D 
        || item_Length!=0);
  
  delete [] buffer;
  if (rflag) return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 

template <class Type>
void DicomHandler::GetValue(char * _val, Type & _rval) {
  
#if BYTE_ORDER == BIG_ENDIAN
  if(fLittleEndian) {      // little endian
#else // BYTE_ORDER == LITTLE_ENDIAN
    if(!fLittleEndian) {     // big endian
#endif
      const G4int SIZE = sizeof(_rval);
      char ctemp;
      for(G4int i = 0; i < SIZE/2; ++i) {
        ctemp = _val[i];
        _val[i] = _val[SIZE - 1 - i];
        _val[SIZE - 1 - i] = ctemp;
      }
    }
    _rval = *(Type *)_val;
  }
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 
