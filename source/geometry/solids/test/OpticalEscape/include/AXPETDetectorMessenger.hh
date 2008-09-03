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
// $Id: AXPETDetectorMessenger.hh,v 1.1 2008-09-03 13:34:03 gcosmo Exp $
// ------------------------------------------------------------
// Geant4 class header file
//
// 03/09/2008, by T.Nikitina
// ------------------------------------------------------------

#ifndef AXPETDetectorMessenger_h
#define AXPETDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class AXPETDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;

class AXPETDetectorMessenger: public G4UImessenger
{
  public:
    AXPETDetectorMessenger(AXPETDetectorConstruction * myDC);
    ~AXPETDetectorMessenger();
      void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    AXPETDetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selDetCmd;
    G4UIcmdWithADouble * rotXCmd;
    G4UIcmdWithADouble * rotYCmd;
    G4UIcmdWithADouble * rotZCmd;
    G4UIcmdWithAnInteger* AbortCmd;
};

#endif


