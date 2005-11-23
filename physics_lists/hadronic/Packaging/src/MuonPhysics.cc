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
 #include "MuonPhysics.hh"
 #include "G4ParticleDefinition.hh"
 #include "G4ParticleTable.hh"

 #include "G4MuonPlus.hh"
 #include "G4MuonMinus.hh"
 #include "G4TauMinus.hh"
 #include "G4TauPlus.hh"
 #include "G4NeutrinoTau.hh"
 #include "G4AntiNeutrinoTau.hh"
 #include "G4NeutrinoMu.hh"
 #include "G4AntiNeutrinoMu.hh"

 #include "G4ProcessManager.hh"

 #include "globals.hh"
 #include "G4ios.hh"

 MuonPhysics::MuonPhysics(const G4String& name)
                    :  G4VPhysicsConstructor(name)
 {
 }

 MuonPhysics::~MuonPhysics()
 {
   if(wasActivated)
   {
     G4ProcessManager * pManager = 0;
//     pManager = G4MuonPlus::MuonPlus()->GetProcessManager();
     pManager = G4MuonMinus::MuonMinus()->GetProcessManager();
     if(pManager) pManager->RemoveProcess(&fMuMinusCaptureAtRest);

   }
 }

 void MuonPhysics::ConstructProcess()
 {
   G4ProcessManager * pManager = 0;

   wasActivated = true;
   // Muon Plus Physics

   // Muon Minus Physics
   pManager = G4MuonMinus::MuonMinus()->GetProcessManager();
    // add capture processes, not in EMstandard
   pManager->AddRestProcess(&fMuMinusCaptureAtRest);

 }

 void MuonPhysics::ConstructParticle()
 {}



 // 2002 by J.P. Wellisch
