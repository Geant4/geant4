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
//---------------------------------------------------------------
//
//  G4NuclideTableMessenger.hh
//
//  Class Description:
//    This is a messenger class to interface to exchange information
//    between ParticleDefinition and UI.
//
//  /particle/manage/nuclide   Nuclide Table control commands.
//   Commands : 
//     lifetime * Set threshold of half-life.
//
//  History:
//    11 November 2015, T. Koi   : The 1st version created.
//
//---------------------------------------------------------------

#ifndef G4NuclideTableMessenger_h
#define G4NuclideTableMessenger_h 1

class G4NuclideTable;

class G4UIdirectory;
//class G4UIcmdWithoutParameter;
//class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
//class G4UIcmdWithAnInteger; 

#include "G4UImessenger.hh"
#include "globals.hh"

class G4NuclideTableMessenger: public G4UImessenger
{
  public:
    G4NuclideTableMessenger(G4NuclideTable* nuclideTable = 0);
    virtual ~G4NuclideTableMessenger();

  public: // With Description
    virtual void SetNewValue(G4UIcommand * command,G4String newValues);

  private:
    G4NuclideTableMessenger(const G4NuclideTableMessenger&):G4UImessenger(){};

  private:
    G4NuclideTable* theNuclideTable;

    G4UIdirectory *             thisDirectory;
    //G4UIcmdWithoutParameter *   dumpCmd;
    G4UIcmdWithADoubleAndUnit * lifetimeCmd; 
    G4UIcmdWithADoubleAndUnit * lToleranceCmd; 
};

#endif
