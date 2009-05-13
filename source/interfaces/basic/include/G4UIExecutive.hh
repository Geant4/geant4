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
// $Id: G4UIExecutive.hh,v 1.1 2009-05-13 09:01:36 kmura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ====================================================================
//   G4UIExecutive.hh
//
//   This class helps automatic instantiation of user session
//   according to your environment variable like G4UI_USE_XXX.
//
//   Usage in main():
//
//   ...
//   #include "G4UIExecutive.hh"
//
//   int main(int argc, char** argv)
//   {
//     ...
//     G4UIExecutive* myapp = new G4UIExecutive(argc, argv);
//     myapp-> SessionStart();
//     ...
//     delete myapp;
//     ...
//
// ====================================================================
#ifndef G4UI_EXECUTIVE_H
#define G4UI_EXECUTIVE_H 1

#include "G4UImanager.hh"
#include "G4UIsession.hh"

#if defined(G4UI_USE_TCSH)
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#elif defined(G4UI_USE_XM)
#include "G4UIXm.hh"

#elif defined(G4UI_USE_WIN32)
#include "G4UIWin32.hh"

#elif defined(G4UI_USE_QT)
#include "G4UIQt.hh"
#include "G4Qt.hh"

#else
#include "G4UIterminal.hh"
#include "G4UIcsh.hh"

#endif

class G4UIExecutive {
private:
  G4UIsession* session;
  G4VUIshell* shell;

  G4String batchname;
  G4bool qbatch;

public:
  G4UIExecutive(G4int argc=1, char** argv=0);
  ~G4UIExecutive();

  void SetPrompt(const G4String& prompt);
  void SetLsColor(TermColorIndex dirColor, TermColorIndex cmdColor);

  void SessionStart();
  
};

#include "G4UIExecutive.icc"

#endif
