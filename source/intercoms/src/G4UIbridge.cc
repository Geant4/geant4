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
// G4UIbridge
//
// Author: A.Dotti, 2013
// --------------------------------------------------------------------
#include "G4UIbridge.hh"
#include "G4UImanager.hh"

// --------------------------------------------------------------------
G4UIbridge::G4UIbridge(G4UImanager* localUI, G4String dir)
  : localUImanager(localUI)
{
  // make sure dirName starts and ends with '/'
  if(dir[0] == '/')
  {
    dirName = dir;
  }
  else
  {
    dirName = "/" + dir;
  }
  if(dirName.back() != '/')
  {
    dirName += "/";
  }

  // register to the master G4UImanager
  G4UImanager* masterUI = G4UImanager::GetMasterUIpointer();
  if(masterUI != nullptr)
  {
    masterUI->RegisterBridge(this);
  }
  else
  {
    G4Exception("G4UIbridge::G4UIbridge()", "UI7001", FatalException,
                "G4UImanager for the master thread is not yet instantiated. "
                "Instantiate G4MTRunManager first.");
  }
}

// --------------------------------------------------------------------
G4int G4UIbridge::ApplyCommand(const G4String& aCmd)
{
  return localUImanager->ApplyCommand(aCmd);
}
