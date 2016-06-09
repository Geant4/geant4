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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4GoudsmitSaundersonTable
//
// Author:        Omrane Kadri
//
// Creation date: 20.02.2009
//
// Modifications:
// 04.03.2009 V.Ivanchenko cleanup and format according to Geant4 EM style
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4GoudsmitSaundersonTable_h
#define G4GoudsmitSaundersonTable_h 1

#include "globals.hh"

class G4GoudsmitSaundersonTable 
{
public:

  G4GoudsmitSaundersonTable();

  ~G4GoudsmitSaundersonTable();

  G4double SampleTheta(G4double , G4double ,G4double);

private:  

  //  hide assignment operator
  G4GoudsmitSaundersonTable & operator=(const  G4GoudsmitSaundersonTable &right);
  G4GoudsmitSaundersonTable(const  G4GoudsmitSaundersonTable&);

  void LoadPDFandCPDFdata();

  static G4double* PDF;
  static G4double* CPDF;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif

