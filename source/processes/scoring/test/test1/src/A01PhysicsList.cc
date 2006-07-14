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
// $Id: A01PhysicsList.cc,v 1.1 2006-07-14 14:43:29 asaim Exp $
// --------------------------------------------------------------
//
// 28-Jan-04 Add QGSP_BERT and QGSP_BIC for hadronic lists. T. Koi
// 22-Nov-04 Comment out QGSP_BERT and QGSP_BIC
//           Output Notificaiton message             
//           All Particles are created in GeneralPhysics 

#include "A01PhysicsList.hh"

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
#include <iomanip>

#include "A01Transport.hh"
#include "A01ParaPhysics.hh"
#include "A01GeneralPhysics.hh"
#include "A01EMPhysics.hh"
#include "A01MuonPhysics.hh"
#include "A01HadronPhysics.hh"
#include "A01IonPhysics.hh"


A01PhysicsList::A01PhysicsList(G4bool ifPara,G4bool ifDisp):  G4VModularPhysicsList()
{
  // default cut value  (1.0mm)
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // Transportation
  RegisterPhysics( new A01Transport("Transport") );

  // Parallel worl navigation process
  if(ifPara) RegisterPhysics( new A01ParaPhysics("parallel navi"));

  // General Physics ( Create ALL Particle and apply Decay )
  RegisterPhysics( new A01GeneralPhysics("general") );

  // EM Physics ( Apply related Processes to gamma and e-/+)
  RegisterPhysics( new A01EMPhysics("standard EM",ifDisp));

  // Muon Physics ( Apply related processes to mu and tau
  RegisterPhysics(  new A01MuonPhysics("muon"));

   // Hadron Physics ( Apply related processes to hadrons )
  RegisterPhysics(  new A01HadronPhysics("hadron"));

  // Ion Physics ( Apply related processes to ions )
  RegisterPhysics( new A01IonPhysics("ion"));

}

A01PhysicsList::~A01PhysicsList()
{
}

void A01PhysicsList::ConstructProcess()
{
  G4PhysConstVector::iterator itr;
  for (itr = physicsVector->begin(); itr!= physicsVector->end(); ++itr) {
    (*itr)->ConstructProcess();
  }
}

void A01PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  SetCutsWithDefault();
}



