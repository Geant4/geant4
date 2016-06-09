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
#ifndef G4RadioactivityTable_h
#define G4RadioactivityTable_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4RadioactivityTable.hh
//
// Version:             0.b.4
// Date:                14/04/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------.
//
// 29 October 2010, F Lei, QinetiQ, UK
// First created     
//
///////////////////////////////////////////////////////////////////////////////

#include "globals.hh"
#include <map>
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"

using std::map;


class G4RadioactivityTable
{
  // class description
  // G4RadioactivityTable is the table of the radioactive isotopes and their
  // activities.  It is used in the variance reduction mode to recorde the 
  // induced radioactivity in the RDM sensitive volumes.
  // class description - end

public:
  // constructor
  //
  G4RadioactivityTable ();

  // destructor
  virtual ~G4RadioactivityTable();

  //  G4RadioactivityTable(const  G4RadioactivityTable &right);
  
public:
  // with description
  //
  void AddIsotope( G4int, G4int, G4double, G4double, G4double);
  // Add an isotope (Z,A,E) of rate to the table.
  //
  G4TwoVector GetRate(G4int, G4int, G4double);
  // Get the rate of isotoe (Z,A,E)
  //
  G4int Entries() const;
  // Get the total number of isotope spieces in the table
  //
  map<G4ThreeVector,G4TwoVector>*  GetTheMap();
  // Return the stl map! Mainly for printing.
  //
private:
  
  map<G4ThreeVector,G4TwoVector> fRadioactivity;
  
};

#endif






