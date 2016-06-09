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
// $Id$
//
// Created by Mikhail Kosov 6-Nov-2009
//
// --------------------------------------------------------------
// Short description: Algorithm of Synchrotron Radiation from PDG
//
// --------------------------------------------------------------

#ifndef G4QSynchRad_h
#define G4QSynchRad_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VDiscreteProcess.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4ThreeVector.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Gamma.hh"


class G4QSynchRad : public G4VDiscreteProcess
{
public:

  G4QSynchRad(const G4String& processName = "CHIPS_SynchrotronRadiation",
                              G4ProcessType type = fElectromagnetic);

  virtual ~G4QSynchRad(){;}

private:

  G4QSynchRad& operator=(const G4QSynchRad &right);

  G4QSynchRad(const G4QSynchRad&);

public:

  G4double GetMeanFreePath(const G4Track& track, G4double step, G4ForceCondition* fCond);
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);
  G4bool IsApplicable(const G4ParticleDefinition& pd) { return (pd.GetPDGCharge() != 0.); }

  void SetMinGamma(G4double ming) {minGamma = ming;}
  G4double GetMinGamma()          {return minGamma;}
  G4double GetRadius(const G4Track& track); // Revolution Radius (independent units)

// Body
private:

  G4double minGamma;
  G4ThreeVector Polarization;  
};


#endif
