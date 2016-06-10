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
// $Id: G4ParticleMessenger.hh 69557 2013-05-08 12:01:40Z gcosmo $
//
//
//---------------------------------------------------------------
//
//  G4ParticleMessenger.hh
//
//  Class Description:
//    This is a messenger class to interface to exchange information
//    between Particle related classes and UI.
//
// ------------------------------------------------------
//  the List of Directory and Commands
// ------------------------------------------------------
// G4ParticleMessenger
//  /particle/   Paricle control commands.
//   Commands : 
//    select * Select particle 
//    list * List name of particles.
//    find * find particle by PDG encoding.
//    verbose * Set Verbose level of Particle Table
//
// G4ParticlePropertyMessenger
//  /particle/property/   Paricle Table control commands.
//   Commands : 
//     dump * dump particle properties.
//     stable * Set stable flag.
//     lifetime * Set life time.
//     verbose * Set Verbose level
//
// G4DecayTableMessenger
//  /particle/property/decay/   Decay Table control commands.
//   Commands : 
//     select * Enter index of decay mode.
//     dump * Dump decay mode information.
//     br * Set branching ratio. [0< BR <1.0]
//
//
//  History:
//    13 June 1997, H. Kurashige   : The 1st version created.
//    10 Nov 1997,  H.Kurashige    : add /particle/property/Verbose 
//    08 Jan. 1998, H. Kurashige   : new UIcommand
//    08 June 1998, H. Kurashige   : remove fProcessManagerMessenger
//    25 Nov. 1998, H. Kurashige   : add /particle/find
//    08 Jun. 2008, H. Kurashige   : add /particle/verbose
//    30 Jul. 2009, H. Kurashige   : add /particle/createAllIon
//    1  May. 2013, H. Kurashige   : add /particle/createAllIsomer
//---------------------------------------------------------------

#ifndef G4ParticleMessenger_h
#define G4ParticleMessenger_h 1

class G4ParticleDefinition;
class G4ParticleTable;
class G4ParticlePropertyMessenger;


class G4UIdirectory;
class G4UIcmdWithAString; 
class G4UIcmdWithAnInteger; 
class G4UIcmdWithoutParameter; 

#include "G4UImessenger.hh"
#include "globals.hh"


class G4ParticleMessenger: public G4UImessenger
{
  public: 
    G4ParticleMessenger(G4ParticleTable* pTable = 0);
    virtual ~G4ParticleMessenger();

  public: // With Description
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
  //  !!!  can not use "copy constructor" !!!!
    G4ParticleMessenger(const G4ParticleMessenger&):G4UImessenger(){};

  private:
    G4UIdirectory *             thisDirectory;
    G4UIcmdWithAString *        listCmd;
    G4UIcmdWithAString *        selectCmd;
    G4UIcmdWithAnInteger *      findCmd;
    G4UIcmdWithoutParameter *   createAllIonCmd;
    G4UIcmdWithoutParameter *   createAllIsomerCmd;
    G4UIcmdWithAnInteger *      verboseCmd;

    G4ParticleTable* theParticleTable;
    G4ParticleDefinition* currentParticle;

    G4ParticlePropertyMessenger* fParticlePropertyMessenger;
};

#endif









