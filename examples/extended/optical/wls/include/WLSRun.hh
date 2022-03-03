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
/// \file optical/wls/include/WLSRun.hh
/// \brief Definition of the WLSRun class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSRun_h
#define WLSRun_h 1

#include "G4Run.hh"

class WLSRun : public G4Run
{
 public:
  WLSRun();
  ~WLSRun();

  void AddTIR(G4int n)
  {
    G4double nd(n);
    fNTIR += nd;
    fNTIR2 += nd * nd;
  };
  void AddExiting(G4int n)
  {
    G4double nd(n);
    fNExiting += nd;
    fNExiting2 += nd * nd;
  };
  void AddEscapedEnd(G4int n)
  {
    G4double nd(n);
    fEscapedEnd += nd;
    fEscapedEnd2 += nd * nd;
  };
  void AddEscapedMid(G4int n)
  {
    G4double nd(n);
    fEscapedMid += nd;
    fEscapedMid2 += nd * nd;
  };
  void AddBounce(G4int n)
  {
    G4double nd(n);
    fBounce += nd;
    fBounce2 += nd * nd;
  };
  void AddWLSBounce(G4int n)
  {
    G4double nd(n);
    fWLSBounce += nd;
    fWLSBounce2 += nd * nd;
  };
  void AddClad1Bounce(G4int n)
  {
    G4double nd(n);
    fClad1Bounce += nd;
    fClad1Bounce2 += nd * nd;
  };
  void AddClad2Bounce(G4int n)
  {
    G4double nd(n);
    fClad2Bounce += nd;
    fClad2Bounce2 += nd * nd;
  };
  void AddReflected(G4int n)
  {
    G4double nd(n);
    fReflected += nd;
    fReflected2 += nd * nd;
  };
  void AddEscaped(G4int n)
  {
    G4double nd(n);
    fEscaped += nd;
    fEscaped2 += nd * nd;
  };
  void AddMirror(G4int n)
  {
    G4double nd(n);
    fMirror += nd;
    fMirror2 += nd * nd;
  };
  void AddDetectorHits(G4int n)
  {
    G4double nd(n);
    fDetectorHits += nd;
    fDetectorHits2 += nd * nd;
  };

  void EndOfRun();
  void Merge(const G4Run*) override;

 private:
  G4double fNTIR;
  G4double fNTIR2;
  G4double fNExiting;
  G4double fNExiting2;
  G4double fEscapedEnd;
  G4double fEscapedEnd2;
  G4double fEscapedMid;
  G4double fEscapedMid2;
  G4double fBounce;
  G4double fBounce2;
  G4double fWLSBounce;
  G4double fWLSBounce2;
  G4double fClad1Bounce;
  G4double fClad1Bounce2;
  G4double fClad2Bounce;
  G4double fClad2Bounce2;
  G4double fReflected;
  G4double fReflected2;
  G4double fEscaped;
  G4double fEscaped2;
  G4double fMirror;
  G4double fMirror2;
  G4double fDetectorHits;
  G4double fDetectorHits2;
};

#endif
