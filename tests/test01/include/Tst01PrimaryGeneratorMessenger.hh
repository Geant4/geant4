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

