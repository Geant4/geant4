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
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsNuBeam
//
// Author: Julia Yarba, FNAL/CD (2014)
// Comment: somewhat "molded" after HadronPhysicsFTFP_BETT
//
// Modified:
// 18.07.2017 A.Dotti: refactoring following new standard
// 02.10.2020 V.Ivanchenko: use more methods from G4HadronPhysicsFTFP_BERT 
//            base class; code clean-up 
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsNuBeam_h
#define G4HadronPhysicsNuBeam_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4HadronPhysicsNuBeam : public G4HadronPhysicsFTFP_BERT
{

public: 
  G4HadronPhysicsNuBeam(G4int verbose =1);
  G4HadronPhysicsNuBeam(const G4String& name, G4bool quasiElastic=false);
  virtual ~G4HadronPhysicsNuBeam() {}

  void ConstructProcess() override;

  // copy constructor and hide assignment operator
  G4HadronPhysicsNuBeam(G4HadronPhysicsNuBeam &) = delete;
  G4HadronPhysicsNuBeam & operator =
  (const G4HadronPhysicsNuBeam &right) = delete;

protected:
  //Modify the minimum needed
  virtual void Proton() override;

private:
  G4double maxFTFP_proton;
};

#endif

