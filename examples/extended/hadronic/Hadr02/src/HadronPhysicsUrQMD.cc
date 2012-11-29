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
/// \file hadronic/Hadr02/src/HadronPhysicsUrQMD.cc
/// \brief Implementation of the HadronPhysicsUrQMD class
//
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2012 A. Dotti
//   created from HadronPhysicsUrQMD
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifdef G4_USE_URQMD
#include "HadronPhysicsUrQMD.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4QHadronInelasticDataSet.hh"
#include "G4ProcessManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronPhysicsUrQMD::HadronPhysicsUrQMD(G4int)
  :  G4VPhysicsConstructor("hInelastic UrQMD")
{
  fNeutrons = 0;
  fUrQMDNeutron = 0;
  fLEPNeutron = 0;
  fPiK = 0;
  fUrQMDPiK = 0;
  fPro = 0;
  fUrQMDPro = 0;    
  fHyperon = 0;
  fAntiBaryon = 0;
  fUrQMDAntiBaryon = 0;
  fCHIPSInelastic = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronPhysicsUrQMD::CreateModels()
{
  fNeutrons=new G4NeutronBuilder;
  fUrQMDNeutron=new UrQMDNeutronBuilder();
  fNeutrons->RegisterMe(fUrQMDNeutron);
  //Need LEP for capture process
  fNeutrons->RegisterMe(fLEPNeutron=new G4LEPNeutronBuilder);
  fLEPNeutron->SetMinInelasticEnergy(0.0*eV);   // no inelastic from LEP
  fLEPNeutron->SetMaxInelasticEnergy(0.0*eV);  

  fPro=new G4ProtonBuilder;
  fUrQMDPro=new UrQMDProtonBuilder();
  fPro->RegisterMe(fUrQMDPro);

  fPiK=new G4PiKBuilder;
  fUrQMDPiK=new UrQMDPiKBuilder();
  fPiK->RegisterMe(fUrQMDPiK);
  
  //For Hyperons use FTF model
  fHyperon=new G4HyperonFTFPBuilder;
    
  fAntiBaryon=new G4AntiBarionBuilder;
  fUrQMDAntiBaryon=new  UrQMDAntiBarionBuilder();
  fAntiBaryon->RegisterMe( fUrQMDAntiBaryon );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronPhysicsUrQMD::~HadronPhysicsUrQMD()
{
  delete fNeutrons;
  delete fUrQMDNeutron;
  delete fLEPNeutron;    

  delete fPiK;
  delete fUrQMDPiK;
    
  delete fPro;
  delete fUrQMDPro;    
    
  delete fHyperon;
  delete fAntiBaryon;
  delete fUrQMDAntiBaryon;
  
  delete fCHIPSInelastic;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronPhysicsUrQMD::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronPhysicsUrQMD::ConstructProcess()
{
  CreateModels();
  fNeutrons->Build();
  fPro->Build();
  fPiK->Build();
  // use CHIPS cross sections also for Kaons
  fCHIPSInelastic = new G4QHadronInelasticDataSet();
  
  FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(fCHIPSInelastic);
  FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(fCHIPSInelastic);
  FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(fCHIPSInelastic);
  FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(fCHIPSInelastic);

  fHyperon->Build();
  fAntiBaryon->Build();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HadronicProcess* 
HadronPhysicsUrQMD::FindInelasticProcess(const G4ParticleDefinition* p)
{
  G4HadronicProcess* had = 0;
  if(p) {
     G4ProcessVector*  pvec = p->GetProcessManager()->GetProcessList();
     size_t n = pvec->size();
     if(0 < n) {
       for(size_t i=0; i<n; ++i) {
         if(fHadronInelastic == ((*pvec)[i])->GetProcessSubType()) {
           had = static_cast<G4HadronicProcess*>((*pvec)[i]);
           break;
         }
       }
     }
  }
  return had;
}

#endif //G4_USE_URQMD
