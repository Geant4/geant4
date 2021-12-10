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
/// \file optical/wls/include/WLSEventAction.hh
/// \brief Definition of the WLSEventAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSEventAction_h
#define WLSEventAction_h 1

#include "G4Types.hh"
#include "G4UserEventAction.hh"

class WLSEventActionMessenger;

class WLSEventAction : public G4UserEventAction
{
 public:
  WLSEventAction();
  ~WLSEventAction();

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  G4int GetEventNo();
  void SetEventVerbose(G4int);

  void AddTIR() { fNTIR += 1; };
  void AddExiting() { fNExiting += 1; };
  void AddEscapedEnd() { fEscapedEnd += 1; };
  void AddEscapedMid() { fEscapedMid += 1; };
  void AddBounce() { fBounce += 1; };
  void AddWLSBounce() { fWLSBounce += 1; };
  void AddClad1Bounce() { fClad1Bounce += 1; };
  void AddClad2Bounce() { fClad2Bounce += 1; };
  void AddReflected() { fReflected += 1; };
  void AddEscaped() { fEscaped += 1; };
  void AddMirror() { fMirror += 1; };

 private:
  WLSEventActionMessenger* fEventMessenger;

  G4int fVerboseLevel;

  G4int fMPPCCollID;

  G4int fNTIR;
  G4int fNExiting;
  G4int fEscapedEnd;
  G4int fEscapedMid;
  G4int fBounce;
  G4int fWLSBounce;
  G4int fClad1Bounce;
  G4int fClad2Bounce;
  G4int fReflected;
  G4int fEscaped;
  G4int fMirror;
};

#endif
