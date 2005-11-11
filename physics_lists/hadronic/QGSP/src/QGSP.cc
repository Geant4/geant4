//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: QGSP.cc,v 1.1 2005-11-11 22:57:16 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   QGSP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard 
//
//----------------------------------------------------------------------------
//

#include "QGSP.hh"
#include "G4DecayBuilder.hh"
#include "G4EmExtraBuilder.hh"
#include "G4EmStandardBuilder.hh"
#include "G4IonBuilder.hh"
#include "G4HadronBuilderQGSP.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

//#include "G4EmPhysicsListMessenger.hh"


QGSP::QGSP(): G4VModularPhysicsList() 
{
  verbose = 0;
  G4LossTableManager::Instance();
  defaultCutValue = 0.7*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  //  pMessenger = new G4EmPhysicsListMessenger(this);

  G4cout << "### QGSP: You are using the simulation engine: QGSP 2.8"<<G4endl;
  G4cout <<G4endl<<G4endl;

  // EM Physics
  this->RegisterPhysics(new G4EmStandardBuilder());
  RegisterPhysics(new G4EmExtraBuilder());

  // General Physics
  RegisterPhysics(new G4DecayBuilder());

  // Hadron Physics
  RegisterPhysics(new G4HadronBuilderQGSP());

  // Ion Physics
  RegisterPhysics(new G4IonBuilder());
}

QGSP::~QGSP()
{
  // delete pMessenger;
}

void QGSP::SetCuts()
{
  if (verbose >1){
    G4cout << "QGSP::SetCuts:" << G4endl;
  }  

  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
  SetParticleCuts(cutForElectron, G4Electron::Electron());
  SetParticleCuts(cutForPositron, G4Positron::Positron());
 
  if(verbose > 0) DumpCutValuesTable();  
}

void QGSP::SetVerbose(G4int val)
{
  verbose = val;
}

void QGSP::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

void QGSP::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

void QGSP::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}
