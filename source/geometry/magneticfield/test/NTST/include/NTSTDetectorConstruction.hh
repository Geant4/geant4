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
// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTDetectorConstruction.hh,v 1.3 2006-06-29 18:25:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef NTSTDetectorConstruction_H
#define NTSTDetectorConstruction_H 1
#include "G4Transform3D.hh"
#include "globals.hh"

class NTSTFileRead;
class G4VPhysicalVolume;
class NTSTDetectorMessenger;
class G4LogicalVolume;
class G4ChordFinder;

#include "G4VUserDetectorConstruction.hh"

#include "NTSTField.hh"


class NTSTDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    NTSTDetectorConstruction();
    ~NTSTDetectorConstruction();
    void SetInputFileName(G4String);
    void SetDebugCmd(G4int);
    void SetOuterRadius(G4double);
    void SetNSubLayer(G4int);
    void PrintCorners(const G4Transform3D&, G4LogicalVolume*);
    void DisableDetector(G4String);
    
public:
    G4VPhysicalVolume* Construct();
    
    void GetFieldCallStats() {
    	G4cout << "Number calls to field = " << field.GetCount() << G4endl;
	field.ClearCount();
    }

private:
    NTSTFileRead* _FileRead;
    G4bool debug;
    G4double radius; // outer radius of the SVT mother volume
    NTSTDetectorMessenger* DetectorMessenger;
    G4int NSubLayer; // default number of layers
    G4bool disableSVT;
    G4bool disableDCH;
    
    NTSTField field;
    G4ChordFinder *fpChordFinder; 
    G4double  fMinChordStep;   // Minimum Step for chord
};

#endif

