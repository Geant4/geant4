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
#include "G4HadronicInteractionWrapper.hh"

G4HadFinalState* G4HadronicInteractionWrapper::
ApplyInteraction(G4HadProjectile& thePro, 
                 G4Nucleus& targetNucleus,
                 G4HadronicInteraction* theInteraction,
                 const G4String& theProcessName,
                 const G4String& theModelName)
{
  static G4ThreadLocal G4HadronicWhiteBoard  *theBoard_G4MT_TLS_ = 0 ; if (!theBoard_G4MT_TLS_) theBoard_G4MT_TLS_ = & G4HadronicWhiteBoard::Instance();G4HadronicWhiteBoard &theBoard = *theBoard_G4MT_TLS_;
  theBoard.SetProjectile(thePro);
  theBoard.SetTargetNucleus(targetNucleus);
  theBoard.SetProcessName(theProcessName);
  theBoard.SetModelName(theModelName);
  return theInteraction->ApplyYourself( thePro, targetNucleus);
}
