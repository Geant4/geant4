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
// $Id: G4ParticlePropertyMessenger.hh 67971 2013-03-13 10:13:24Z gcosmo $
//
//
//---------------------------------------------------------------
//
//  G4ParticlePropertyMessenger.hh
//
//  Class Description:
//    This is a messenger class to interface to exchange information
//    between ParticleDefinition and UI.
//
//  /particle/property/   Paricle Table control commands.
//   Commands : 
//     dump * dump particle properties.
//     stable * Set stable flag.
//     lifetime * Set life time.
//     Verbose * Set Verbose level
//
//  History:
//    13 June 1997, H. Kurashige   : The 1st version created.
//    13 Nov. 1997, H. Kurashige   : fix bugs
//    08 Jan. 1998, H. Kurashige   : new UIcommand
//    08 Apr. 1999, H. Kurashige   : fix some improper codings
//
//---------------------------------------------------------------

#ifndef G4ParticlePropertyMessenger_h
#define G4ParticlePropertyMessenger_h 1

class G4ParticleTable;
class G4ParticleDefinition;
class G4DecayTableMessenger;

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger; 

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ParticlePropertyMessenger: public G4UImessenger
{
  public:
    G4ParticlePropertyMessenger(G4ParticleTable* pTable = 0);
    virtual ~G4ParticlePropertyMessenger();

  public: // With Description
    virtual void SetNewValue(G4UIcommand * command,G4String newValues);
    virtual G4String GetCurrentValue(G4UIcommand * command);

  private:
    G4ParticlePropertyMessenger(const G4ParticlePropertyMessenger&):G4UImessenger(){};

  private:
    G4ParticleDefinition* SetCurrentParticle();
    G4ParticleTable* theParticleTable;
    G4ParticleDefinition* currentParticle;

    G4UIdirectory *             thisDirectory;
    G4UIcmdWithoutParameter *   dumpCmd;
    G4UIcmdWithABool *          stableCmd; 
    G4UIcmdWithAnInteger *      verboseCmd;
    G4UIcmdWithADoubleAndUnit * lifetimeCmd; 
 
    G4DecayTableMessenger* fDecayTableMessenger;
};

#endif






