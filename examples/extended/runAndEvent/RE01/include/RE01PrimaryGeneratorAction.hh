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
/// \file runAndEvent/RE01/include/RE01PrimaryGeneratorAction.hh
/// \brief Definition of the RE01PrimaryGeneratorAction class
//
//

#ifndef RE01PrimaryGeneratorAction_h
#define RE01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4VPrimaryGenerator;
class G4Event;
class RE01PrimaryGeneratorMessenger;

class RE01PrimaryGeneratorAction:public G4VUserPrimaryGeneratorAction
{
public:
  RE01PrimaryGeneratorAction();
  virtual ~RE01PrimaryGeneratorAction();

public:
  virtual void GeneratePrimaries(G4Event* anEvent);
  inline  void SetHEPEvtGenerator(G4bool f)
  { fUseHEPEvt = f; }
  inline G4bool GetHEPEvtGenerator()
  { return fUseHEPEvt; }

private:
  G4VPrimaryGenerator* fHEPEvt;
  G4VPrimaryGenerator* fParticleGun;
  RE01PrimaryGeneratorMessenger* fMessenger;
  G4bool fUseHEPEvt;

};

#endif


