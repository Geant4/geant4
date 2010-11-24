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
// $Id: Tst65DetectorMessenger.hh,v 1.1 2010-11-24 15:15:29 tkoi Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst65DetectorMessenger_h
#define Tst65DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst65DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst65DetectorMessenger: public G4UImessenger
{
  public:
    Tst65DetectorMessenger(Tst65DetectorConstruction * myDC);
   ~Tst65DetectorMessenger();
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst65DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selMatCmd;

};

#endif

