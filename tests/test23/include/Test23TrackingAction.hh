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
// $Id: Test23TrackingAction.hh,v 1.2 2006-06-29 21:49:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- Test23TrackingAction header ----------------
//                 by Mikhail Kossov, December 2003.
//  Test23TrackingAction class of the CHIPS Test of G4QCaptureAtRest process in GEANT4
// -----------------------------------------------------------------------------------
// It collects secondaries and measures the mean energy of the each sort of particles
// -----------------------------------------------------------------------------------

#ifndef Test23TrackingAction_h
#define Test23TrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTypes.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

class Test23TrackingAction : public G4UserTrackingAction
{

public:
  Test23TrackingAction();
  ~Test23TrackingAction();

  void PreUserTrackingAction(const G4Track*);
  void PostUserTrackingAction(const G4Track*) {}; // inline void void

  void PrintResult();
  void ResetResult();

private:

  Test23TrackingAction & operator=(const Test23TrackingAction &right);
  Test23TrackingAction(const Test23TrackingAction&);

private:

  static G4int    n_Gammas;
  static G4int    n_NeutMu;
  static G4int    n_AnNuEl;
  static G4int    nElectrs;
  static G4int    nPositrs;
  static G4int    nProtons;
  static G4int    nNeutron;
  static G4int    n_Deutrs;
  static G4int    n_Heliu3;
  static G4int    n_Triton;
  static G4int    n_Alphas;
  static G4int    n_Others;
  static G4double E_Gammas;
  static G4double E_NeutMu;
  static G4double E_AnNuEl;
  static G4double EElectrs;
  static G4double EPositrs;
  static G4double EProtons;
  static G4double ENeutron;
  static G4double E_Deutrs;
  static G4double E_Heliu3;
  static G4double E_Triton;
  static G4double E_Alphas;
  static G4double E_Others;
};

#endif

