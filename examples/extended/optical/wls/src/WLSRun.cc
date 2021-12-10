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
/// \file optical/wls/src/WLSRun.cc
/// \brief Implementation of the WLSRun class
//
//

#include "WLSRun.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSRun::WLSRun()
  : G4Run()
{
  fNTIR          = 0.;
  fNTIR2         = 0.;
  fNExiting      = 0.;
  fNExiting2     = 0.;
  fEscapedEnd    = 0.;
  fEscapedEnd2   = 0.;
  fEscapedMid    = 0.;
  fEscapedMid2   = 0.;
  fBounce        = 0.;
  fBounce2       = 0.;
  fWLSBounce     = 0.;
  fWLSBounce2    = 0.;
  fClad1Bounce   = 0.;
  fClad1Bounce2  = 0.;
  fClad2Bounce   = 0.;
  fClad2Bounce2  = 0.;
  fReflected     = 0.;
  fReflected2    = 0.;
  fEscaped       = 0.;
  fEscaped2      = 0.;
  fMirror        = 0.;
  fMirror2       = 0.;
  fDetectorHits  = 0.;
  fDetectorHits2 = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSRun::~WLSRun() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSRun::Merge(const G4Run* run)
{
  const WLSRun* localRun = static_cast<const WLSRun*>(run);

  fNTIR += localRun->fNTIR;
  fNTIR2 += localRun->fNTIR2;
  fNExiting += localRun->fNExiting;
  fNExiting2 += localRun->fNExiting2;
  fEscapedEnd += localRun->fEscapedEnd;
  fEscapedEnd2 += localRun->fEscapedEnd2;
  fEscapedMid += localRun->fEscapedMid;
  fEscapedMid2 += localRun->fEscapedMid2;
  fBounce += localRun->fBounce;
  fBounce2 += localRun->fBounce2;
  fWLSBounce += localRun->fWLSBounce;
  fWLSBounce2 += localRun->fWLSBounce2;
  fClad1Bounce += localRun->fClad1Bounce;
  fClad1Bounce2 += localRun->fClad1Bounce2;
  fClad2Bounce += localRun->fClad2Bounce;
  fClad2Bounce2 += localRun->fClad2Bounce2;
  fReflected += localRun->fReflected;
  fReflected2 += localRun->fReflected2;
  fEscaped += localRun->fEscaped;
  fEscaped2 += localRun->fEscaped2;
  fMirror += localRun->fMirror;
  fMirror2 += localRun->fMirror2;
  fDetectorHits += localRun->fDetectorHits;
  fDetectorHits2 += localRun->fDetectorHits2;

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSRun::EndOfRun()
{
  if(numberOfEvent == 0)
    return;
  G4double TotNbofEvents = G4double(numberOfEvent);

  fNTIR           = fNTIR / TotNbofEvents;
  fNTIR2          = fNTIR2 / TotNbofEvents;
  G4double rmsTIR = fNTIR2 - fNTIR * fNTIR;
  if(rmsTIR > 0.)
    rmsTIR = std::sqrt(rmsTIR);
  else
    rmsTIR = 0.;

  fNExiting           = fNExiting / TotNbofEvents;
  fNExiting2          = fNExiting2 / TotNbofEvents;
  G4double rmsExiting = fNExiting2 - fNExiting * fNExiting;
  if(rmsExiting > 0.)
    rmsExiting = std::sqrt(rmsExiting);
  else
    rmsExiting = 0.;

  fEscapedEnd            = fEscapedEnd / TotNbofEvents;
  fEscapedEnd2           = fEscapedEnd2 / TotNbofEvents;
  G4double rmsEscapedEnd = fEscapedEnd2 - fEscapedEnd * fEscapedEnd;
  if(rmsEscapedEnd > 0.)
    rmsEscapedEnd = std::sqrt(rmsEscapedEnd);
  else
    rmsEscapedEnd = 0.;

  fEscapedMid            = fEscapedMid / TotNbofEvents;
  fEscapedMid2           = fEscapedMid2 / TotNbofEvents;
  G4double rmsEscapedMid = fEscapedMid2 - fEscapedMid * fEscapedMid;
  if(rmsEscapedMid > 0.)
    rmsEscapedMid = std::sqrt(rmsEscapedMid);
  else
    rmsEscapedMid = 0.;

  fBounce            = fBounce / TotNbofEvents;
  fBounce2           = fBounce2 / TotNbofEvents;
  G4double rmsBounce = fBounce2 - fBounce * fBounce;
  if(rmsBounce > 0.)
    rmsBounce = std::sqrt(rmsBounce);
  else
    rmsBounce = 0.;

  fWLSBounce            = fWLSBounce / TotNbofEvents;
  fWLSBounce2           = fWLSBounce2 / TotNbofEvents;
  G4double rmsWLSBounce = fWLSBounce2 - fWLSBounce * fWLSBounce;
  if(rmsWLSBounce > 0.)
    rmsWLSBounce = std::sqrt(rmsWLSBounce);
  else
    rmsWLSBounce = 0.;

  fClad1Bounce            = fClad1Bounce / TotNbofEvents;
  fClad1Bounce2           = fClad1Bounce2 / TotNbofEvents;
  G4double rmsClad1Bounce = fClad1Bounce2 - fClad1Bounce * fClad1Bounce;
  if(rmsClad1Bounce > 0.)
    rmsClad1Bounce = std::sqrt(rmsClad1Bounce);
  else
    rmsClad1Bounce = 0.;

  fClad2Bounce            = fClad2Bounce / TotNbofEvents;
  fClad2Bounce2           = fClad2Bounce2 / TotNbofEvents;
  G4double rmsClad2Bounce = fClad2Bounce2 - fClad2Bounce * fClad2Bounce;
  if(rmsClad2Bounce > 0.)
    rmsClad2Bounce = std::sqrt(rmsClad2Bounce);
  else
    rmsClad2Bounce = 0.;

  fReflected            = fReflected / TotNbofEvents;
  fReflected2           = fReflected2 / TotNbofEvents;
  G4double rmsReflected = fReflected2 - fReflected * fReflected;
  if(rmsReflected > 0.)
    rmsReflected = std::sqrt(rmsReflected);
  else
    rmsReflected = 0.;

  fEscaped            = fEscaped / TotNbofEvents;
  fEscaped2           = fEscaped2 / TotNbofEvents;
  G4double rmsEscaped = fEscaped2 - fEscaped * fEscaped;
  if(rmsEscaped > 0.)
    rmsEscaped = std::sqrt(rmsEscaped);
  else
    rmsEscaped = 0.;

  fMirror            = fMirror / TotNbofEvents;
  fMirror2           = fMirror2 / TotNbofEvents;
  G4double rmsMirror = fMirror2 - fMirror * fMirror;
  if(rmsMirror > 0.)
    rmsMirror = std::sqrt(rmsMirror);
  else
    rmsMirror = 0.;

  fDetectorHits            = fDetectorHits / TotNbofEvents;
  fDetectorHits2           = fDetectorHits2 / TotNbofEvents;
  G4double rmsDetectorHits = fDetectorHits2 - fDetectorHits * fDetectorHits;
  if(rmsDetectorHits > 0.)
    rmsDetectorHits = std::sqrt(rmsDetectorHits);
  else
    rmsDetectorHits = 0.;

  G4int prec = G4cout.precision(3);
  G4cout << "\n ======================== run summary =====================\n";
  G4cout << "Average number per event:" << G4endl;
  G4cout << " TIR:          " << fNTIR << " +- " << rmsTIR << G4endl
         << " Exiting:      " << fNExiting << " +- " << rmsExiting << G4endl
         << " Escaped Mid:  " << fEscapedMid << " +- " << rmsEscapedMid
         << G4endl << " Escaped End:  " << fEscapedEnd << " +- "
         << rmsEscapedEnd << G4endl << " Bounced:      " << fBounce << " +- "
         << rmsBounce << G4endl << " WLS Bounce:   " << fWLSBounce << " +- "
         << rmsWLSBounce << G4endl << " Clad1 Bounce: " << fClad1Bounce
         << " +- " << rmsClad1Bounce << G4endl
         << " Clad2 Bounce: " << fClad2Bounce << " +- " << rmsClad2Bounce
         << G4endl << " Reflected:    " << fReflected << " +- " << rmsReflected
         << G4endl << " Escaped:      " << fEscaped << " +- " << rmsEscaped
         << G4endl << " Mirror:       " << fMirror << " +- " << rmsMirror
         << G4endl << " Detector hit: " << fDetectorHits << " +- "
         << rmsDetectorHits << G4endl;

  G4cout << G4endl;
  G4cout.precision(prec);
}
