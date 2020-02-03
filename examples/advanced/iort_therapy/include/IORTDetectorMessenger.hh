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
//   (e) University of Wollongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#ifndef IORTDetectorMessenger_h
#define IORTDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class IORTDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithoutParameter; 
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;        

class IORTDetectorMessenger: public G4UImessenger
{
public:
  IORTDetectorMessenger(IORTDetectorConstruction* );
  ~IORTDetectorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:

  // Pointer to the phantom/detector 
  IORTDetectorConstruction* iortDetector;

  G4UIdirectory *changeThePhantomDir,  *changeTheDetectorDir, *changeTheDisc1Dir, *changeTheDisc2Dir, *deleteTheDiscDir, *insertTheDiscDir, *changeTheAnglediscDir;

  G4UIcmdWithoutParameter   *updateCmd, *deletediscCmd, *insertdiscCmd;
  G4UIcmdWithAString        *changeThePhantomMaterialCmd, *changeTheDisc1MaterialCmd, *changeTheDisc2MaterialCmd; 
  G4UIcmdWith3VectorAndUnit *changeThePhantomSizeCmd,
    *changeThePhantomPositionCmd, 
    *changeTheDetectorSizeCmd, 
    *changeTheDetectorToPhantomPositionCmd,
    *changeTheDetectorVoxelCmd;

  G4UIcmdWithADoubleAndUnit   *changeOuterRadiusDiscoIORTCmd, *changeinnerRadiusDiscoIORTCmd, *changeheightDiscoIORTCmd,*changeDiscoXPositionIORTCmd, *changeDiscoYPositionIORTCmd, *changeDiscoZPositionIORTCmd, *changeOuterRadiusDisco1IORTCmd, *changeinnerRadiusDisco1IORTCmd, *changeheightDisco1IORTCmd,*changeDisco1XPositionIORTCmd, *changeTheAnglediscCmd;
};
#endif

