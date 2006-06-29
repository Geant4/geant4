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

#ifndef Tst01PrimaryGeneratorMessenger_h
#define Tst01PrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst01PrimaryGeneratorAction ;
class G4UIcmdWithoutParameter ;
class G4UIcmdWithAString ;
class G4UIcmdWith3VectorAndUnit ;

class Tst01PrimaryGeneratorMessenger: public G4UImessenger
{
  public:

    Tst01PrimaryGeneratorMessenger(Tst01PrimaryGeneratorAction * myPGA);

    void SetNewValue(G4UIcommand * command,G4String newValues);

  private:

    Tst01PrimaryGeneratorAction * fPrimaryGeneratorAction;

    G4UIcmdWithoutParameter * standardGun ;
    G4UIcmdWithAString * randomGun ;
    G4UIcmdWith3VectorAndUnit* viewerCmd ;
    G4UIcmdWith3VectorAndUnit* planeCmd ;
};

#endif

