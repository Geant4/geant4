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
//
// DicomConfiguration.hh :
// - Handling of the header of *.g4 files
// - Reading <Data.dat> file

#ifndef DicomConfiguration_h
#define DicomConfiguration_h 1

#include "globals.hh"
#include <vector>

class DicomConfiguration
{
public:

   DicomConfiguration();
  ~DicomConfiguration() { }

  // This function reads <Data.dat>, return true or false
  // if successfull or not
  G4bool ReadDataFile();

    G4int GetCompressionValue() {return compressionValue;}

  G4int GetTotalNumberOfFile() {return totalNumberOfFile;}

  std::vector<G4String> GetListOfFile() {return listOfFile;}

  G4int GetTotalRows() {return totalRows;}

  G4int GetTotalColumns() {return totalColumns;}  

    G4int GetTotalPixels() {return totalPixels;}  
    //G4int GetTotalPixels() {return densityValue.size();}

  G4double GetXPixelSpacing() {return  xPixelSpacing;}

  G4double GetYPixelSpacing() {return  yPixelSpacing;}

  G4double GetSliceThickness() {return sliceTickness;} 

  std::vector<G4double> GetSliceLocation() {return  sliceLocation;} 

  G4int IsCompressionUsed() {return compressionUsed;}

  G4double GetDensityValue(G4int i);

    void ClearDensityData() {densityValue.clear(); };

    G4int ReadG4File(G4String g4File);
 
private:
    
    static short compressionValue;
    static G4int totalNumberOfFile;
    static std::vector<G4String> listOfFile;
    static short totalRows;
    static short totalColumns;
    static G4int totalPixels;
    static G4double xPixelSpacing;
    static G4double yPixelSpacing;
    static G4double sliceTickness;
    static std::vector<G4double> sliceLocation; 
    static short compressionUsed;
    static std::vector<G4double> densityValue;	

};
#endif

