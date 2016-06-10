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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy
// (This class was adapted from the purging magnet example, by S.Larsson and J. Generowicz.)
#include "globals.hh"
#include "G4ElectricField.hh"
#include "G4ElectroMagneticField.hh"
#include "G4ios.hh"

#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

class HadrontherapyElectricTabulatedField3D

 : public G4ElectricField

{
  
  // Storage space for the table
  vector< vector< vector< G4double > > > xEField;
  vector< vector< vector< G4double > > > yEField;
  vector< vector< vector< G4double > > > zEField;
  // The dimensions of the table
  G4int Enx,Eny,Enz; 
  // The physical limits of the defined region
  G4double Eminx, Emaxx, Eminy, Emaxy, Eminz, Emaxz;
  // The physical extent of the defined region
  G4double dx1, dy1, dz1;
  G4double feXoffset;
  G4double feYoffset;
  G4double feZoffset;
  G4bool einvertX, einvertY, einvertZ;

public:
  HadrontherapyElectricTabulatedField3D(const char* filename, G4double exOffset, G4double eyOffset, G4double ezOffset );
  void  GetFieldValue( const  G4double Epoint[4],
		       G4double *Efield) const;
};

