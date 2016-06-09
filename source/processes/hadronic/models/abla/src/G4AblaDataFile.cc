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
// $Id$ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#include "G4AblaDataFile.hh"
//#include "G4HadronicException.hh"
#include "globals.hh" // Needed for G4Exception.
#include <fstream>

G4AblaDataFile::G4AblaDataFile()
{
  verboseLevel = 0;
}

/**
 * Read all data from files.
 */
bool G4AblaDataFile::readData()
{
  if(!getenv("G4ABLADATA")) {
    //    throw G4HadronicException(__FILE__, __LINE__, "ERROR: Data
    //    missing. Set environment variable G4ABLA3.0 to point to the
    //    directory containing data files needed by INCL and ABLA
    //    models.");
    G4String errorMessage1 = "ERROR: Data missing. Set environment variable G4ABLADATA\n";
    G4String errorMessage2 = "\t to point to the directory containing data files needed\n";
    G4String errorMessage3 = "\t by INCL and ABLA models.\n";
    G4String errorMessage = errorMessage1 + errorMessage2 + errorMessage3;
    G4Exception(errorMessage);
  }
  
  G4String dataPath(getenv("G4ABLADATA"));
  G4String flAlphaFile(dataPath + "/flalpha.dat");
  G4String frldmFile(  dataPath + "/frldm.dat");
  G4String vgsldFile(  dataPath + "/vgsld.dat");
  G4String pace2File(  dataPath + "/pace2.dat");

  if(verboseLevel > 1) {
    G4cout <<"Data path   = " << dataPath    << G4endl;
    G4cout <<"FlAlphaFile = " << flAlphaFile << G4endl;
    G4cout <<"FrldmFile   = " << frldmFile   << G4endl;
    G4cout <<"VgsldFile   = " << vgsldFile   << G4endl;
    G4cout <<"Pace2File   = " << pace2File   << G4endl;
  }
  
  std::ifstream flalphain(flAlphaFile.c_str());
  std::ifstream frldmin(frldmFile.c_str());  
  std::ifstream vgsldin(vgsldFile.c_str());  
  std::ifstream pace2in(pace2File.c_str());

  std::filebuf *buf1 = flalphain.rdbuf();
  std::filebuf *buf2 = frldmin.rdbuf();
  std::filebuf *buf3 = vgsldin.rdbuf();
  std::filebuf *buf4 = pace2in.rdbuf();  
  if (!((buf1->is_open()) && (buf2->is_open()) && (buf3->is_open()) && (buf4->is_open()))) {
    G4Exception("ERROR: Data missing. Could not find ABLA data file in " + dataPath +
		" defined by environment variable G4ABLADATA");
  }
  
  G4double flalpha, frldm, vgsld, pace2;
  const G4int rows = 98;
  const G4int cols = 153;
  const G4int massnumbers = 263;
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      setAlpha(j, i, 0.0);
      setEcnz( j, i, 0.0);
      setVgsld(j, i, 0.0);      
    }
  }
  
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      flalphain >> flalpha;
      frldmin >> frldm;
      vgsldin >> vgsld;      
      setAlpha(j, i, flalpha);
      setEcnz( j, i, frldm);
      setVgsld(j, i, vgsld);      
    }
  }
  flalphain.close();
  frldmin.close();  
  vgsldin.close();

  int A = 0, Zbegin = 0, Zend = 0;
  G4String str1, str2, str3;
  for(int i = 0; i < 500; i++) {
    for(int j = 0; j < 500; j++) {
      setPace2(i, j, 0.0);
    }
  }
  
  for(int i = 0; i < massnumbers; i++) {
    pace2in >> str1 >> A >> str2 >> Zbegin >> str3 >> Zend;
    for(int j = Zbegin; j <= Zend; j++) {
      pace2in >> pace2;
      setPace2(A, j, pace2);
    }
  }
  pace2in.close();
  if(std::fabs(getPace2(A, Zend) - 114516.10) > 1e-6) {
    G4cout <<"ERROR: Problem in parsing datafile " + pace2File << G4endl;
    return false;
  }

  return true;
}


