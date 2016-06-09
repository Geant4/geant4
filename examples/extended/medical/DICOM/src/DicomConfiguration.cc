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
// DicomConfiguration.cc :
//	- Handling of the header of *.g4 files
// 	- Reading <Data.dat> file

#include "globals.hh"
#include "DicomConfiguration.hh"
#include <fstream>
#include <vector>

short DicomConfiguration::compressionValue = 0;
G4int DicomConfiguration::totalNumberOfFile = 0;
std::vector<G4String> DicomConfiguration::listOfFile;
short DicomConfiguration::totalRows = 0;
short DicomConfiguration::totalColumns = 0;
G4int DicomConfiguration::totalPixels = 0;
G4double DicomConfiguration::xPixelSpacing = 0.;
G4double DicomConfiguration::yPixelSpacing = 0.;
G4double DicomConfiguration::sliceTickness = 0.;
std::vector<G4double> DicomConfiguration::sliceLocation;
short DicomConfiguration::compressionUsed = 0;
std::vector<G4double> DicomConfiguration::densityValue;

//
DicomConfiguration::DicomConfiguration() {
    ReadDataFile();
}

//
G4bool DicomConfiguration::ReadDataFile() {

    if(totalNumberOfFile > 0) return true;

    totalPixels = 0;

    std::ifstream dataFile("Data.dat");
    G4String nameOfFileBuffer;
    if(dataFile.good() != 1 ) return 1;

    dataFile >> compressionValue;
    dataFile >> totalNumberOfFile;
 
    for(G4int i = 0; i < totalNumberOfFile; i++ ) {

	dataFile >> nameOfFileBuffer;
	listOfFile.push_back( nameOfFileBuffer );

	// read densities from .g4 file
	ReadG4File(nameOfFileBuffer);
    }

    dataFile.close();
    return 0;
}

G4int DicomConfiguration::ReadG4File( G4String g4File ) {

    //densityValue.clear();

    g4File = g4File + ".g4";
    std::ifstream readingG4FileHeader(g4File.c_str(),
				      std::ios_base::in | std::ios_base::binary);

    if ( readingG4FileHeader.good() != 1 ) return 1;

    readingG4FileHeader.read((char *)&totalRows, 2);
    readingG4FileHeader.read((char *)&totalColumns, 2);
    readingG4FileHeader.read((char *)&xPixelSpacing, 8);
    readingG4FileHeader.read((char *)&yPixelSpacing, 8);
    readingG4FileHeader.read((char *)&sliceTickness, 8);
    G4double sliceLocationBuff;
    readingG4FileHeader.read((char *)&sliceLocationBuff, 8);
    readingG4FileHeader.read((char *)&compressionUsed, 2);

    sliceLocation.push_back(sliceLocationBuff);

    G4double density;
    for(int y = 0; y < totalRows/compressionUsed; y++) {
	for(int x = 0; x < totalColumns/compressionUsed; x++) {
	    readingG4FileHeader.read((char *)&density, sizeof(G4double));
	    densityValue.push_back(density);
	    totalPixels++;
	}
    }

    readingG4FileHeader.close();
    return 0;
}

G4double DicomConfiguration::GetDensityValue(G4int i) {

    G4double value = 0.;

    if (i >= 0) {
	unsigned int j = i;
	//
	if(j > densityValue.size() ) {
	    // Throw exception, return dummy, cerr error message...
	    G4cout << "out of range in GetDensityValue()! : "
		   << j << ", " << totalPixels << G4endl;
	} else {
	    value = densityValue[i];
	}
    } else {
	G4cout << "out of range in GetDensityValue()!"<<G4endl;
    }

    return value;
}
