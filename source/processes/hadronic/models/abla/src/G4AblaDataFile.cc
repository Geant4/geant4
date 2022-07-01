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
// ABLAXX statistical de-excitation model
// Jose Luis Rodriguez, GSI (translation from ABLA07 and contact person)
// Pekka Kaitaniemi, HIP (initial translation of ablav3p)
// Aleksandra Kelic, GSI (ABLA07 code)
// Davide Mancusi, CEA (contact person INCL)
// Aatos Heikkinen, HIP (project coordination)
//

#define ABLAXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4AblaDataFile.hh"

#ifdef ABLAXX_IN_GEANT4_MODE
#include "globals.hh"
#else
#include "G4INCLGeant4Compat.hh"
#endif

#include <fstream>
#include <cmath>
#include <iostream>
#include <cstdlib>

#ifdef ABLAXX_IN_GEANT4_MODE
G4AblaDataFile::G4AblaDataFile() {
#else
G4AblaDataFile::G4AblaDataFile(G4INCL::Config *config)
  : G4AblaVirtualData(config) {
  theConfig = config;
#endif
  verboseLevel = 0;
}

G4AblaDataFile::~G4AblaDataFile()
{
}

/**
 * Read all data from files.
 */
bool G4AblaDataFile::readData()
{
#ifdef ABLAXX_IN_GEANT4_MODE
  if(!G4FindDataDir("G4ABLADATA")) {
    //    throw G4HadronicException(__FILE__, __LINE__, "ERROR: Data
    //    missing. Set environment variable G4ABLA3.0 to point to the
    //    directory containing data files needed by INCL and ABLA
    //    models.");
    // G4String errorMessage1 = "ERROR: Data missing. Set environment variable G4ABLADATA\n";
    // G4String errorMessage2 = "\t to point to the directory containing data files needed\n";
    // G4String errorMessage3 = "\t by INCL and ABLA models.\n";
    // G4String errorMessage = errorMessage1 + errorMessage2 + errorMessage3;
    // G4Exception(errorMessage);
    G4ExceptionDescription ed;
    ed << " Data missing: set environment variable G4ABLADATA\n"
       << " to point to the directory containing data files needed\n"
       << " by the ABLA model" << G4endl;
    G4Exception("G4AblaDataFile::readData()","ABLA_001",
                FatalException, ed);
  }
  G4String dataPath(G4FindDataDir("G4ABLADATA"));
#else
  G4String dataPath(theConfig->getABLAXXDataFilePath().c_str());
#endif
  G4String flAlphaFile(dataPath + "/flalpha.dat");
  G4String frldmFile(  dataPath + "/frldm.dat");
  G4String vgsldFile(  dataPath + "/vgsld.dat");
  G4String pace2File(  dataPath + "/pace2.dat");
  G4String rmsFile(  dataPath + "/rms.dat");
  G4String defoFile(  dataPath + "/defo.dat");
  G4String massFile(  dataPath + "/mass2003.dat");

  if(verboseLevel > 1) {
    // G4cout <<"Data path   = " << dataPath    << G4endl;
    // G4cout <<"FlAlphaFile = " << flAlphaFile << G4endl;
    // G4cout <<"FrldmFile   = " << frldmFile   << G4endl;
    // G4cout <<"VgsldFile   = " << vgsldFile   << G4endl;
    // G4cout <<"Pace2File   = " << pace2File   << G4endl;
  }
  
  std::ifstream flalphain(flAlphaFile.c_str());
  std::ifstream frldmin(frldmFile.c_str());  
  std::ifstream vgsldin(vgsldFile.c_str());  
  std::ifstream pace2in(pace2File.c_str());
  std::ifstream rmsin(rmsFile.c_str());
  std::ifstream defoin(defoFile.c_str());
  std::ifstream massin(massFile.c_str());

  std::filebuf *buf1 = flalphain.rdbuf();
  std::filebuf *buf2 = frldmin.rdbuf();
  std::filebuf *buf3 = vgsldin.rdbuf();
  std::filebuf *buf4 = pace2in.rdbuf();  
  std::filebuf *buf5 = rmsin.rdbuf(); 
  std::filebuf *buf6 = defoin.rdbuf();
  std::filebuf *buf7 = massin.rdbuf();
  if (!((buf1->is_open()) && (buf2->is_open()) && (buf3->is_open()) && (buf4->is_open()) && (buf5->is_open()) && (buf6->is_open()) && (buf7->is_open()))) {
#ifdef ABLAXX_IN_GEANT4_MODE
    G4ExceptionDescription ed;
    ed << "Data missing: could not find ABLA data file in " << dataPath
       << "defined by environment variable G4ABLADATA" << G4endl;
    G4Exception("G4AblaDataFile::readData()", "ABLA", FatalException, ed);
#else
    std::cerr << "Error opening file." << std::endl;
#endif
  }
  
  G4double fflalpha, ffrldm, fvgsld, fpace2, frms;
  int fj,fk,a2,a3,a4;
  G4double fbeta2,fbeta4;
  G4double a7;
  const G4int rows = 99;
  const G4int cols = 154;
  const G4int rowsbeta = 137;
  const G4int colsbeta = 251;
  const G4int rowsmass = 13;
  const G4int colsmass = 154;
  const G4int massnumbers = 263;
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      setAlpha(j, i, 0.0);
      setEcnz( j, i, 0.0);
      setVgsld(j, i, 0.0);
      setRms(j, i, 0.0);
    }
  }
  
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      flalphain >> fflalpha;
      frldmin >> ffrldm;
      vgsldin >> fvgsld;  
      rmsin >> frms;    
      setAlpha(j, i, fflalpha);
      setEcnz( j, i, ffrldm);
      setVgsld(j, i, fvgsld);
      setRms(j, i, frms);
    }
  }


  for(int i = 0; i < rowsbeta; i++) {
    for(int j = 0; j < colsbeta; j++) {
      setBeta2(j, i, 0.0);
      setBeta4(j, i, 0.0);
    }
  }

  for(int i = 0; i < 8983; i++) {
  defoin >> fj >> fk >> fbeta2 >> fbeta4;
  setBeta2(fk, fj, fbeta2);
  setBeta4(fk, fj, fbeta4);
  }

  for(int i = 0; i < rowsmass; i++) {
    for(int j = 0; j < colsmass; j++) {
      setMexp(j, i, 0.0);
      setMexpID(j,i,0);
    }
  }
   massin >> a2 >> a3 >> a4 >> a7 ;
   while(!massin.eof()){
   //
     if(a3<13.){
     setMexpID(a2,a3,1);
     setMexp(a2,a3,938.7829835*a3+939.5653301*a2-1.*a4*a7/1000.);
     }
     massin >> a2 >> a3 >> a4 >> a7 ;
   }

  flalphain.close();
  frldmin.close();  
  vgsldin.close();
  rmsin.close();
  defoin.close();
  massin.close();

  G4String str1, str2, str3;
  for(int i = 0; i < 500; i++) {
    for(int j = 0; j < 500; j++) {
      setPace2(i, j, 0.0);
    }
  }
  
  int A = 0, Zbegin = 0, Zend = 0;
  for(int i = 0; i < massnumbers; i++) {
    pace2in >> str1 >> A >> str2 >> Zbegin >> str3 >> Zend;
    if(Zbegin >= 0 && Zbegin < getPaceCols() &&
       A >= 0 && A < getPaceRows()) {
      for(int j = Zbegin; j <= Zend; j++) {
	pace2in >> fpace2;
	setPace2(A, j, fpace2);
      }
    }
  }
  pace2in.close();
  if(std::abs(getPace2(A, Zend) - 114516.10) > 1e-6) {
    std::cerr << "ERROR: Problem in parsing datafile " + pace2File << std::endl;
    return false;
  }

  return true;
}

