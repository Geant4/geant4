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
//
// $Id: DetectorConstruction.hh,v 1.1 2008/11/20 15:41:54 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// Class DetectorConstruction
//
// A detector construction class loading the geometry from GDML files.
//
// ----------------------------------------------------------------------------

#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1


#include "G4VUserDetectorConstruction.hh"

#include "G4Material.hh"

#include "globals.hh"

#include "G4GDMLParser.hh"

class DetectorMessenger;

// ----------------------------------------------------------------------------

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    // Constructor and destructor
    //
    DetectorConstruction();
   ~DetectorConstruction();

    // Construction of Detector
    //
    G4VPhysicalVolume* Construct();

    // Make List of materials
    //
    void ListOfMaterials();

    // Reading GDML
    //
    void SetReadFile( const G4String& File );

  private:

    G4Material* Air ;
    G4Material* Aluminum ;
    G4Material* Pb;
    G4Material* Xenon;

    // Extended reader
    //
    G4GDMLReadStructure* reader;

    // GDMLparser
    //
    G4GDMLParser* parser;
        
    // Reading Settings
    //
    G4String fReadFile;

    // Detector Messenger
    //
    DetectorMessenger* detectorMessenger;
};

// ----------------------------------------------------------------------------

#endif
