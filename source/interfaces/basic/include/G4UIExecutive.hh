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
// $Id: G4UIExecutive.hh 66892 2013-01-17 10:57:59Z gunter $
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
//     if (session->IsGUI())
//       // Do any extra for a GUI session
//
//     myapp-> SessionStart();
//     ...
//     delete myapp;
//     ...
//
// ====================================================================
#ifndef G4UI_EXECUTIVE_HH
#define G4UI_EXECUTIVE_HH

#include "G4VUIshell.hh"
#include <map>

class G4UIsession;

class G4UIExecutive {
public:
  G4UIExecutive(G4int argc, char** argv, const G4String& type = "");
  ~G4UIExecutive();

  G4UIsession* GetSession() const;

  G4bool IsGUI() const;

  void SetPrompt(const G4String& prompt);
  void SetLsColor(TermColorIndex dirColor, TermColorIndex cmdColor);

  void SessionStart();

private:
  enum SessionType { kNone, kQt, kXm, kWin32, kWt, kGag, kTcsh, kCsh };
  SessionType selected;

  G4UIsession* session;
  G4VUIshell* shell;

  G4bool isGUI;

  std::map<G4String, G4String> sessionMap;

  void SelectSessionByArg(const G4String& stype);
  void SelectSessionByEnv();
  void SelectSessionByFile(const G4String& appname);
  void SelectSessionByBestGuess();
};

// ====================================================================
inline G4UIsession* G4UIExecutive::GetSession() const
{
  return session;
}

inline G4bool G4UIExecutive::IsGUI() const
{
  return isGUI;
}

#endif
