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
// G4UIsession
//
// Author: Makoto Asai, 1998
// --------------------------------------------------------------------

#include "G4UIsession.hh"

G4int G4UIsession::inSession = 0;

// --------------------------------------------------------------------
G4UIsession::G4UIsession()
{
  ++inSession;
}

// --------------------------------------------------------------------
G4UIsession::G4UIsession(G4int iBatch)
  : ifBatch(iBatch)
{
}

// --------------------------------------------------------------------
G4UIsession::~G4UIsession()
{
  if(ifBatch == 0)
  {
    --inSession;
  }
}

// --------------------------------------------------------------------
G4UIsession* G4UIsession::SessionStart()
{
  return nullptr;
}

// --------------------------------------------------------------------
void G4UIsession::PauseSessionStart(const G4String&)
{
}

// --------------------------------------------------------------------
G4int G4UIsession::InSession()
{
  return inSession;
}

// --------------------------------------------------------------------
G4int G4UIsession::ReceiveG4cout(const G4String& coutString)
{
  std::cout << coutString << std::flush;
  return 0;
}

// --------------------------------------------------------------------
G4int G4UIsession::ReceiveG4cerr(const G4String& cerrString)
{
  std::cerr << cerrString << std::flush;
  return 0;
}
