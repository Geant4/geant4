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
// $Id: G4AblaVirtualData.hh,v 1.1 2008-02-27 18:31:11 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#ifndef G4AblaVirtualData_hh
#define G4AblaVirtualData_hh 1

#include "globals.hh"

/**
 * An interface to data used by INCL and ABLA. This interface allows
 * us to abstract the actual source of data. Currently the data is
 * read from datafiles by using class G4InclAblaDataFile.  @see
 * G4InclAblaDataFile
 */

class G4AblaVirtualData {
protected:

  /**
   * Constructor
   */
  G4AblaVirtualData();

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

  G4double getAlpha(G4int A, G4int Z);

  /**
   * Get the value of Alpha.
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

  G4int getAlphaRows();
  G4int getAlphaCols();

  G4int getPaceRows();
  G4int getPaceCols();

  virtual G4bool readData() = 0;
	
private:

  static const G4int alphaRows = 155;
  static const G4int alphaCols = 100;

  static const G4int paceRows = 500;
  static const G4int paceCols = 500;

  G4double alpha[alphaRows][alphaCols];
  G4double ecnz[alphaRows][alphaCols];
  G4double vgsld[alphaRows][alphaCols];
  G4double pace2[paceRows][paceCols];
};

#endif
