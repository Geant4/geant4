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
// Jose Luis Rodriguez, UDC (translation from ABLA07 and contact person)
// Pekka Kaitaniemi, HIP (initial translation of ablav3p)
// Aleksandra Kelic, GSI (ABLA07 code)
// Davide Mancusi, CEA (contact person INCL)
// Aatos Heikkinen, HIP (project coordination)
//

#include "G4AblaDataFile.hh"

#include "G4AblaDataDefs.hh"
#include "globals.hh"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

G4AblaDataFile::G4AblaDataFile()
{
  verboseLevel = 0;
}

/**
 * Read all data from files.
 */
G4bool G4AblaDataFile::readData()
{
  if (!G4FindDataDir("G4ABLADATA")) {
    G4ExceptionDescription ed;
    ed << " Data missing: set environment variable G4ABLADATA\n"
       << " to point to the directory containing data files needed\n"
       << " by the ABLA model" << G4endl;
    G4Exception("G4AblaDataFile::readData()", "ABLA_001", FatalException, ed);
  }
  G4String dataPath(G4FindDataDir("G4ABLADATA"));

  G4String flAlphaFile(dataPath + "/flalpha.dat");
  G4String frldmFile(dataPath + "/frldm.dat");
  G4String vgsldFile(dataPath + "/vgsld.dat");
  G4String rmsFile(dataPath + "/rms.dat");
  G4String defoFile(dataPath + "/defo.dat");
  G4String massFile(dataPath + "/mass2020.dat");

  if (verboseLevel > 1) {
    // G4cout <<"Data path   = " << dataPath    << G4endl;
    // G4cout <<"FlAlphaFile = " << flAlphaFile << G4endl;
    // G4cout <<"FrldmFile   = " << frldmFile   << G4endl;
    // G4cout <<"VgsldFile   = " << vgsldFile   << G4endl;
  }

  std::ifstream flalphain(flAlphaFile.c_str());
  std::ifstream frldmin(frldmFile.c_str());
  std::ifstream vgsldin(vgsldFile.c_str());
  std::ifstream rmsin(rmsFile.c_str());
  std::ifstream defoin(defoFile.c_str());
  std::ifstream massin(massFile.c_str());

  if (!massin.is_open()) {
    massFile = dataPath + "/mass2016.dat";
    massin.close();
    massin.open(massFile.c_str());
    std::cout << "Mass evaluation file mass2020.dat not found, current file: " << massFile.c_str()
              << std::endl;

    if (!massin.is_open()) {
      massFile = dataPath + "/mass2003.dat";
      massin.close();
      massin.open(massFile.c_str());
      std::cout << "Mass evaluation file mass2016.dat not found, current file: " << massFile.c_str()
                << std::endl;
    }
  }

  std::filebuf* buf1 = flalphain.rdbuf();
  std::filebuf* buf2 = frldmin.rdbuf();
  std::filebuf* buf3 = vgsldin.rdbuf();
  std::filebuf* buf4 = rmsin.rdbuf();
  std::filebuf* buf5 = defoin.rdbuf();
  std::filebuf* buf6 = massin.rdbuf();
  if (!((buf1->is_open()) && (buf2->is_open()) && (buf3->is_open()) && (buf4->is_open())
        && (buf5->is_open()) && (buf6->is_open())))
  {
    G4ExceptionDescription ed;
    ed << "Data missing: could not find ABLA data file in " << dataPath
       << "defined by environment variable G4ABLADATA" << G4endl;
    G4Exception("G4AblaDataFile::readData()", "ABLA", FatalException, ed);
  }

  G4double fflalpha, ffrldm, fvgsld, frms;
  G4int fj = 0, fk = 0, a2, a3, a4;
  G4double fbeta2, fbeta4;
  G4double a7;
  const G4int rows = 99;
  const G4int cols = 154;
  const G4int rowsbeta = 137;
  const G4int colsbeta = 251;

  for (G4int i = 0; i < zcols; i++) {
    for (G4int j = 0; j < nrows; j++) {
      setAlpha(j, i, 0.0);
      setEcnz(j, i, 0.0);
      setVgsld(j, i, 0.0);
      setRms(j, i, 0.0);
    }
  }

  for (G4int i = 0; i < rows; i++) {
    for (G4int j = 0; j < cols; j++) {
      flalphain >> fflalpha;
      frldmin >> ffrldm;
      vgsldin >> fvgsld;
      rmsin >> frms;
      setAlpha(j, i, fflalpha);
      setEcnz(j, i, ffrldm);
      setVgsld(j, i, fvgsld);
      setRms(j, i, frms);
    }
  }

  for (G4int i = 0; i < rowsbeta; i++) {
    for (G4int j = 0; j < colsbeta; j++) {
      setBeta2(j, i, 0.0);
      setBeta4(j, i, 0.0);
    }
  }

  defoin >> fj >> fk >> fbeta2 >> fbeta4;
  while (!defoin.eof()) {
    setBeta2(fk, fj, fbeta2);
    setBeta4(fk, fj, fbeta4);
    defoin >> fj >> fk >> fbeta2 >> fbeta4;
  }

  for (G4int i = 0; i < zcols; i++) {
    for (G4int j = 0; j < nrows; j++) {
      setMexp(j, i, 0.0);
      setMexpID(j, i, 0);
    }
  }
  massin >> a2 >> a3 >> a4 >> a7;
  while (!massin.eof()) {
    //
    if (a3 < lpcols) {
      setMexpID(a2, a3, 1);
      setMexp(a2, a3, 938.7829835 * a3 + 939.5653301 * a2 - 1. * a4 * a7 / 1000.);
    }
    massin >> a2 >> a3 >> a4 >> a7;
  }

  flalphain.close();
  frldmin.close();
  vgsldin.close();
  rmsin.close();
  defoin.close();
  massin.close();

  return true;
}
