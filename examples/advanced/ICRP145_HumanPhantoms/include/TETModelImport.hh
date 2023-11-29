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
// Author: Haegin Han
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Geant4 Contributors: J. Allison and S. Guatelli
//
#ifndef TETModelImport_h
#define TETModelImport_h 1

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>

#include "G4UIExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "G4Tet.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Colour.hh"

// *********************************************************************
// This class imports the phantom data from *.ele, *.node, and
// *.material files.
// -- DataRead: Construct G4Tet by reading data from *.ele and *.node
//              files
// -- MaterialRead: Construct G4Material by reading material data from
//                  *.material file
// -- ColourRead: Construct std::map that contains G4Colour data
//                according to data in colour.dat file
// -- PrintMaterialInformation: Print a table that contains organ ID,
//                              number of tetrahedrons, volume, density,
//                              mass, and name for each organ
// *********************************************************************

class TETModelImport
{
public:
 TETModelImport(G4bool isAF, G4UIExecutive* ui);
 virtual ~TETModelImport() {};
 
 // get methods
 G4String      GetPhantomName()           { return fPhantomName; };
 G4Material*   GetMaterial(G4int idx)     { return fMaterialMap[idx];}
 G4int         GetNumTetrahedron()        { return fTetVector.size();}
 G4int         GetMaterialIndex(G4int idx){ return fMaterialVector[idx]; }
 G4Tet*        GetTetrahedron(G4int idx)  { return fTetVector[idx]; }
 G4double      GetVolume(G4int idx)       { return fVolumeMap[idx]; }
 std::map<G4int, G4double> GetMassMap()   { return fMassMap; }
 std::map<G4int, G4Colour> GetColourMap() { return fColourMap; }
 G4ThreeVector GetPhantomSize()           { return fPhantomSize; }
 G4ThreeVector GetPhantomBoxMin()         { return fBoundingBox_Min; }
 G4ThreeVector GetPhantomBoxMax()         { return fBoundingBox_Max; }

private:
 // private methods
 void DataRead(G4String, G4String);
 void MaterialRead(G4String);
 void ColourRead();
 void PrintMaterialInfomation();

 G4String fPhantomDataPath;
 G4String fPhantomName;

 G4ThreeVector fBoundingBox_Min;
 G4ThreeVector fBoundingBox_Max;
 G4ThreeVector fPhantomSize;

 std::vector<G4ThreeVector> fVertexVector;
 std::vector<G4Tet*>        fTetVector;
 std::vector<G4int*>        fEleVector;
 std::vector<G4int>         fMaterialVector;
 std::map<G4int, G4int>     fNumTetMap;
 std::map<G4int, G4double>  fVolumeMap;
 std::map<G4int, G4double>  fMassMap;
 std::map<G4int, G4Colour>  fColourMap;

 std::map<G4int, std::vector<std::pair<G4int, G4double>>> fMaterialIndexMap;
 std::vector<G4int>                                       fMaterialIndex;
 std::map<G4int, G4Material*>                             fMaterialMap;
 std::map<G4int, G4double>                                fDensityMap;
 std::map<G4int, G4String>                                fOrganNameMap;
};

#endif
