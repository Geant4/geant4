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
// + Universit� Laval, Qu�bec (QC) Canada
//*******************************************************//

//*******************************************************
//
// DicomHandler.hh :
//	- Handling of DICM images
//	- Transforming *.dcm to *.g4 ( pixels->density )
//	- Reading headers and pixels
//	- Transforming pixel to density and creating *.g4
//	  files
//	- Functions are in DicomHandler.cc
//
//
// Base on previous code by :
//	Dragan Tubic <tdragan@gel.ulaval.ca>
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


#ifndef DicomHandler_h
#define DicomHandler_h 1

#include <cstdio>
#include <map>
#include <fstream>

#include "globals.hh"


class DicomHandler
{
public:

  DicomHandler();

    ~DicomHandler();

  G4int ReadFile(FILE *,char *);
  G4int ReadData(FILE *,char *); // note: always use readHeader 
                                    // before readData

  // use ImageMagick to display the image
  //G4int displayImage(char[500]);

    void CheckFileFormat(G4String directory, G4String fileName, G4String calibrationDensityFileName);

private:
    template <class Type> void GetValue(char *, Type &);

private:

    const int DATABUFFSIZE;
    const int LINEBUFFSIZE;
    const int FILENAMESIZE;

    void GetInformation(G4int &, char *);
    G4float Pixel2density(G4int pixel);
    void ReadMaterialIndices( std::ifstream& finData);
    unsigned int GetMaterialIndex( G4float density );
    void StoreData(std::ofstream& foutG4DCM);

    short compression;
    G4int nFiles;
    short rows;
    short columns;
    short bitAllocated;
    G4int maxPixelValue, minPixelValue;
    
    G4double pixelSpacingX, pixelSpacingY;
    G4double sliceThickness;
    G4double sliceLocation;
    
    G4int rescaleIntercept, rescaleSlope;
    
    G4bool littleEndian, implicitEndian;
    short pixelRepresentation;

    G4int** tab;
    std::map<G4float,G4String> fMaterialIndices;

	G4String dicomDirectory, fileData, calibrationDensityFileName;

};
#endif

