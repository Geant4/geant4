//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTDetectorConstruction.hh,v 1.2 2003-12-09 15:35:11 gunter Exp $
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

