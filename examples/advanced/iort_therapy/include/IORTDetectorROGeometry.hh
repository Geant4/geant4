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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wallongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////
//The detectior is devided in voxels. 
//
//
#ifndef IORTDetectorROGeometry_h
#define IORTDetectorROGeometry_h 

#include "G4VReadOutGeometry.hh"

class IORTDetectorROGeometry : public G4VReadOutGeometry
{
public:
  IORTDetectorROGeometry(G4String aString,
				  G4ThreeVector detectorPos,
				  G4double detectorDimX,
				  G4double detectorDimY,
				  G4double detectorDimZ,
				  G4int numberOfVoxelsX,
				  G4int numberOfVoxelsY,
				  G4int numberOfVoxelsZ);

  ~IORTDetectorROGeometry();

private:
  G4VPhysicalVolume* Build();

private:  
  const G4ThreeVector detectorToWorldPosition; 
  const G4double detectorSizeX;
  const G4double detectorSizeY; 
  const G4double detectorSizeZ;

  const G4int numberOfVoxelsAlongX;
  const G4int numberOfVoxelsAlongY; 
  const G4int numberOfVoxelsAlongZ; 
  
  G4VPhysicalVolume *RODetectorZDivisionPhys;
};
#endif
