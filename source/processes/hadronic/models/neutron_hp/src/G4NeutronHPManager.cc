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
// Class Description
// Manager of NetronHP
// 
// 121031 First implementation done by T. Koi (SLAC/PPA)

#include "G4NeutronHPManager.hh"

G4NeutronHPManager* G4NeutronHPManager::instance = NULL;

G4NeutronHPManager::G4NeutronHPManager()
:RWB(NULL)
{
;
}
G4NeutronHPManager::~G4NeutronHPManager()
{
;
}
void G4NeutronHPManager::OpenReactionWhiteBoard()
{
   if ( RWB != NULL ) {
      G4cout << "Warning: G4NeutronHPReactionWhiteBoard is tried doubly opening" << G4endl;
      RWB = new G4NeutronHPReactionWhiteBoard();
   }
   
   RWB = new G4NeutronHPReactionWhiteBoard();
}
G4NeutronHPReactionWhiteBoard* G4NeutronHPManager::GetReactionWhiteBoard()
{
   if ( RWB == NULL ) {
      G4cout << "Warning: try to access G4NeutronHPReactionWhiteBoard before opening" << G4endl;
      RWB = new G4NeutronHPReactionWhiteBoard();
   }
   return RWB; 
}
