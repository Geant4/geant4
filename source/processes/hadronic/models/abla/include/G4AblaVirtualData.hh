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
// Jose Luis Rodriguez, CEA (translation from ABLA07 and contact person)
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Davide Mancusi, CEA (contact person INCL)
// Aatos Heikkinen, HIP (project coordination)
//
#define ABLAXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4AblaVirtualData_hh
#define G4AblaVirtualData_hh 1

#ifdef ABLAXX_IN_GEANT4_MODE
#include "globals.hh"
#else
#include "G4INCLGeant4Compat.hh"
#include "G4INCLConfig.hh"
#endif


/**
 * An interface to data used by ABLA. This interface allows
 * us to abstract the actual source of data. Currently the data is
 * read from datafiles by using class G4AblaDataFile.  @see
 * G4AblaDataFile
 */

class G4AblaVirtualData {
protected:

  /**
   * Constructor, destructor
   */
#ifdef ABLAXX_IN_GEANT4_MODE
  G4AblaVirtualData();
#else
  G4AblaVirtualData(G4INCL::Config *);
#endif
  virtual ~G4AblaVirtualData();

public:
  /**
   * Set the value of Alpha.
   */
  G4bool setAlpha(G4int A, G4int Z, G4double value);

  /**
   * Set the value of Ecnz.
   */
  G4bool setEcnz(G4int A, G4int Z, G4double value);

  /**
   * Set the value of Vgsld.
   */
  G4bool setVgsld(G4int A, G4int Z, G4double value);

  /**
   * Set the value of Pace2.
   */
  G4bool setPace2(G4int A, G4int Z, G4double value);

  /**
   * Set the value of RMS.
   */
  G4bool setRms(G4int A, G4int Z, G4double value);

  /**
   * Set the value of experimental masses.
   */
  G4bool setMexp(G4int A, G4int Z, G4double value);

  /**
   * Set the value of experimental masses ID.
   */
  G4bool setMexpID(G4int A, G4int Z, G4int value);

  /**
   * Set the value of beta2 deformation.
   */
  G4bool setBeta2(G4int A, G4int Z, G4double value);

  /**
   * Set the value of beta4 deformation.
   */
  G4bool setBeta4(G4int A, G4int Z, G4double value);


  /**
   * Get the value of Alpha.
   */
  G4double getAlpha(G4int A, G4int Z);

  /**
   * Get the value of Ecnz.
   */
  G4double getEcnz(G4int A, G4int Z);

  /**
   * Get the value of Vgsld.
   */
  G4double getVgsld(G4int A, G4int Z);

  /**
   * Get the value of Pace2.
   */
  G4double getPace2(G4int A, G4int Z);

  /**
   * Get the value of RMS.
   */
  G4double getRms(G4int A, G4int Z);

  /**
   * Get the value of experimental masses.
   */
  G4double getMexp(G4int A, G4int Z);

  /**
   * Get the value of experimental masses ID.
   */
  G4int getMexpID(G4int A, G4int Z);

  /**
   * Get the value of beta2 deformation.
   */
  G4double getBeta2(G4int A, G4int Z);

  /**
   * Get the value of beta4 deformation.
   */
  G4double getBeta4(G4int A, G4int Z);

  G4int getAlphaRows();
  G4int getAlphaCols();

  G4int getPaceRows();
  G4int getPaceCols();

  virtual G4bool readData() = 0;
	
private:

  static const G4int alphaRows = 154;
  static const G4int alphaCols = 99;

  static const G4int paceRows = 500;
  static const G4int paceCols = 500;

  static const G4int rmsRows = 154;
  static const G4int rmsCols = 99;

  static const G4int betaRows = 251;
  static const G4int betaCols = 137;

  static const G4int massRows = 154;
  static const G4int massCols = 13;

  G4double alpha[alphaRows][alphaCols];
  G4double ecnz[alphaRows][alphaCols];
  G4double vgsld[alphaRows][alphaCols];
  G4double pace2[paceRows][paceCols];
  G4double rms[rmsRows][rmsCols];
  G4double mexp[massRows][massCols];
  G4int mexpid[massRows][massCols];
  G4double beta2[betaRows][betaCols];
  G4double beta4[betaRows][betaCols];
};

#endif
