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
// $Id: G4UserPhysicsListMessenger.hh 66892 2013-01-17 10:57:59Z gunter $
//
// 
//---------------------------------------------------------------
//
//  G4UserPhysicsListMessenger.hh
//
//  Class Description:
//    This is a messenger class to interface to exchange information
//    between ParticleUserList and UI.
// --
//  the List of Directory and Commands
// -       
//  /run/particle/   Paricle control commands.
//   Commands : 
//    SetCuts *  Set default cut value
//    dumpList * Dump List of particles in G4VUserPhysicsList.
//    verbose * Set the Verbose level of G4VUserPhysicsList.
//    addProcessManager * add process manager
//    buildPhysicsTable * build physics table
//    storePhysicsTable * store physics table into files
//    retreivePhysicsTable * retreive physics table from files
//    setStoredInAscii * Switch on/off ascii mode in store/retreive Physics Table
// ------------------------------------------------------------
//	History
//        first version                   09 Jan. 1998 by H.Kurashige 
//        second version                  24 Jan. 1998 by H.Kurashige 
//        add buildPhysicsTable command   13 Apr. 1999 by H.Kurashige
//        add store/retreivePhysicsTable  08 Nov. 2000 by H.Kurashige
//        add setStoredInAscii command    12 Mar. 2001 by H.Kurashige
//        add applyCuts command            2 Aug. 2001 by H.Kurashige
//        add dumpOrderingParam command    3 May. 2011 by H.Kurashige
//        add getCutForAGivenParticle     11 June 2011 by H.Kurashige
// ------------------------------------------------------------

#ifndef G4UserPhysicsListMessenger_h
#define G4UserPhysicsListMessenger_h 1

class G4VUserPhysicsList;

class G4VUserPhysicsList;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString; 
class G4UIcommand;

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UserPhysicsListMessenger: public G4UImessenger
{
  private:
  // hide default constructor
    G4UserPhysicsListMessenger(){}

  public:
    G4UserPhysicsListMessenger(G4VUserPhysicsList* pParticleList);
    virtual ~G4UserPhysicsListMessenger();
    
public: // with description
    virtual  void SetNewValue(G4UIcommand * command,G4String newValues);
    virtual  G4String GetCurrentValue(G4UIcommand * command);

  protected:
    G4VUserPhysicsList* thePhysicsList;
    
  private: //commands
    G4UIdirectory *             theDirectory;
    G4UIcmdWithADoubleAndUnit * setCutCmd; 
    G4UIcommand *               setCutRCmd;
    G4UIcommand *               setCutForAGivenParticleCmd;
    G4UIcmdWithAString *        getCutForAGivenParticleCmd;
    G4UIcmdWithAnInteger *      verboseCmd;
    G4UIcmdWithoutParameter *   dumpListCmd;
    G4UIcmdWithAString *        addProcManCmd;
    G4UIcmdWithAString *        buildPTCmd;
    G4UIcmdWithAString *        storeCmd;
    G4UIcmdWithAString *        retrieveCmd;
    G4UIcmdWithAnInteger *      asciiCmd;
    G4UIcommand *               applyCutsCmd;
    G4UIcmdWithAString *        dumpCutValuesCmd;
    G4UIcmdWithAnInteger*       dumpOrdParamCmd;
};

#endif


