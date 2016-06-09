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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// ParticleSourceMessenger header
// --------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
// This particle source is a shortened version of G4GeneralParticleSource by
// C Ferguson, F Lei & P Truscott (University of Southampton / DERA), with
// some minor modifications.
//////////////////////////////////////////////////////////////////////////////

#ifndef DMXParticleSourceMessenger_h
#define DMXParticleSourceMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DMXParticleSource;

class G4ParticleTable;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;


class DMXParticleSourceMessenger: public G4UImessenger {
  
   public:
     DMXParticleSourceMessenger(DMXParticleSource *fPtclGun);
     ~DMXParticleSourceMessenger();
  
     void SetNewValue(G4UIcommand *command, G4String newValues);
 
  
   private:
     DMXParticleSource *fParticleGun;
     G4ParticleTable *particleTable;
    
   private:
     G4UIdirectory              *gunDirectory;

     G4UIcmdWithAString         *typeCmd;
     G4UIcmdWithAString         *shapeCmd;
     G4UIcmdWith3VectorAndUnit  *centreCmd;
     G4UIcmdWithADoubleAndUnit  *halfzCmd;
     G4UIcmdWithADoubleAndUnit  *radiusCmd;
     G4UIcmdWithAString         *confineCmd;         
     G4UIcmdWithAString         *angtypeCmd;
     G4UIcmdWithAString         *energytypeCmd;
     G4UIcmdWithAnInteger       *verbosityCmd;
     G4UIcommand                *ionCmd;
     G4UIcmdWithAString         *particleCmd;
     G4UIcmdWith3VectorAndUnit  *positionCmd;
     G4UIcmdWith3Vector         *directionCmd;
     G4UIcmdWithADoubleAndUnit  *energyCmd;
     G4UIcmdWithoutParameter    *listCmd;


   private:
     G4bool   fShootIon; 
     G4int    fAtomicNumber;
     G4int    fAtomicMass;
     G4int    fIonCharge;
     G4double fIonExciteEnergy;
  
};

#endif
