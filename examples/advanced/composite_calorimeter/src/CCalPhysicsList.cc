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
//---------------------------------------------------
// This file should be copied as CCalPhysicsList.cc
// if you want to use the 
//                        *** QGSC *** 
// Physics List
//---------------------------------------------------
#include "CCalPhysicsList.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   

#include "GeneralPhysics.hh"
#include "EMPhysics.hh"
#include "MuonPhysics.hh"
#include "HadronPhysicsQGSC.hh"
#include "IonPhysics.hh"

CCalPhysicsList::CCalPhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  // defaultCutValue = 1.0*mm;
  defaultCutValue = 0.7*mm;  
  SetVerboseLevel(1);

  // General Physics
  RegisterPhysics( new GeneralPhysics("general") );

  // EM Physics
  RegisterPhysics( new EMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics(  new MuonPhysics("muon"));

   // Hadron Physics
  RegisterPhysics(  new HadronPhysicsQGSC("hadron"));

  // Ion Physics
  RegisterPhysics( new IonPhysics("ion"));


}

CCalPhysicsList::~CCalPhysicsList()
{
}

void CCalPhysicsList::SetCuts()
{
  if (verboseLevel >1){
    G4cout << "CCalPhysicsList::SetCuts:";
  }  
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 

  SetCutsWithDefault();   
  
  SetCutValue(kNuCut,"nu_e");
  SetCutValue(kNuCut,"nu_tau");
  SetCutValue(kNuCut,"nu_mu");
  SetCutValue(kNuCut,"anti_nu_e");
  SetCutValue(kNuCut,"anti_nu_tau");
  SetCutValue(kNuCut,"anti_nu_mu");
 
  DumpCutValues(G4Electron::ElectronDefinition() ); 
  DumpCutValues(G4Positron::PositronDefinition()); 
  DumpCutValues(G4MuonPlus::MuonPlusDefinition() ); 
  DumpCutValues(G4MuonMinus::MuonMinusDefinition()); 
  DumpCutValues(G4Gamma::GammaDefinition() ) ;
 
  
}



// 2002 by J.P. Wellisch
