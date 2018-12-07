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
/// \file persistency/gdml/G03/include/G03DetectorConstruction.hh
/// \brief Definition of the G03DetectorConstruction class
//
//
//
// Class G03DetectorConstruction
//
// A detector construction class loading the geometry from GDML files.
//
// ----------------------------------------------------------------------------

#ifndef G03DetectorConstruction_H
#define G03DetectorConstruction_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4GDMLParser.hh"

class G4Material;
class G03DetectorMessenger;

/// Detector construction for the GDML extensions example

class G03DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    // Constructor and destructor
    //
    G03DetectorConstruction();
   ~G03DetectorConstruction();

    // Construction of Detector
    //
    virtual G4VPhysicalVolume* Construct();

    // Make List of materials
    //
    void ListOfMaterials();

    // Reading/writing GDML
    //
    void SetReadFile( const G4String& fname );
    void SetWriteFile( const G4String& fname );

  private:

    G4Material* fAir ;
    G4Material* fAluminum ;
    G4Material* fPb;
    G4Material* fXenon;

    // Extended reader
    //
    G4GDMLReadStructure* fReader;

    // Extended writer
    //
    G4GDMLWriteStructure* fWriter;

    // GDMLparser
    //
    G4GDMLParser* fParser;
        
    // Read/write Settings
    //
    G4String fReadFile, fWriteFile;
    G4bool fWritingChoice;
 
    // Detector Messenger
    //
    G03DetectorMessenger* fDetectorMessenger;
};

// ----------------------------------------------------------------------------

#endif
