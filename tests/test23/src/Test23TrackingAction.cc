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
// $Id: Test23TrackingAction.cc,v 1.2 2006-06-29 21:50:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- Test23TrackingAction class ----------------
//                 by Mikhail Kossov, December 2003.
//  Test23TrackingAction class of the CHIPS Simulation Branch in GEANT4

//#define debug

#include "Test23TrackingAction.hh"

G4int Test23TrackingAction::n_Gammas=0;
G4int Test23TrackingAction::n_NeutMu=0;
G4int Test23TrackingAction::n_AnNuEl=0;
G4int Test23TrackingAction::nElectrs=0;
G4int Test23TrackingAction::nPositrs=0;
G4int Test23TrackingAction::nProtons=0;
G4int Test23TrackingAction::nNeutron=0;
G4int Test23TrackingAction::n_Deutrs=0;
G4int Test23TrackingAction::n_Heliu3=0;
G4int Test23TrackingAction::n_Triton=0;
G4int Test23TrackingAction::n_Alphas=0;
G4int Test23TrackingAction::n_Others=0;
G4double Test23TrackingAction::E_Gammas=0;
G4double Test23TrackingAction::E_AnNuEl=0;
G4double Test23TrackingAction::E_NeutMu=0;
G4double Test23TrackingAction::EElectrs=0;
G4double Test23TrackingAction::EPositrs=0;
G4double Test23TrackingAction::EProtons=0;
G4double Test23TrackingAction::ENeutron=0;
G4double Test23TrackingAction::E_Deutrs=0;
G4double Test23TrackingAction::E_Heliu3=0;
G4double Test23TrackingAction::E_Triton=0;
G4double Test23TrackingAction::E_Alphas=0;
G4double Test23TrackingAction::E_Others=0;

Test23TrackingAction::Test23TrackingAction() 
{
#ifdef debug
  G4cout<<"Test23TrackingAction::Constructor is called"<<G4endl;
  PrintResult();
#endif
}

Test23TrackingAction::~Test23TrackingAction() {}

void Test23TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  static G4ParticleDefinition* gamma    = G4Gamma::Gamma();
  static G4ParticleDefinition* neutrino = G4NeutrinoMu::NeutrinoMu();
  static G4ParticleDefinition* antineut = G4AntiNeutrinoE::AntiNeutrinoE();
  static G4ParticleDefinition* electron = G4Electron::Electron();
  static G4ParticleDefinition* positron = G4Positron::Positron();
  static G4ParticleDefinition* proton   = G4Proton::Proton();
  static G4ParticleDefinition* neutron  = G4Neutron::Neutron();
  static G4ParticleDefinition* deuteron = G4Deuteron::Deuteron();
  static G4ParticleDefinition* helium3  = G4He3::He3();
  static G4ParticleDefinition* triton   = G4Triton::Triton();
  static G4ParticleDefinition* alpha    = G4Alpha::Alpha();
  const G4DynamicParticle* secondary = track->GetDynamicParticle();
  const G4ParticleDefinition* particle=secondary->GetDefinition();
  const G4double energy=secondary->GetKineticEnergy();
#ifdef debug
  G4cout<<"Test23TrackingAction::PreUserTrackAct:PDG="<<particle->GetPDGEncoding()
        <<", name="<<particle->GetParticleName()<<G4endl;
#endif
  if     (particle==gamma)    {n_Gammas++; E_Gammas+=energy;}
  else if(particle==neutrino) {n_NeutMu++; E_NeutMu+=energy;}
  else if(particle==antineut) {n_AnNuEl++; E_AnNuEl+=energy;}
  else if(particle==electron) {nElectrs++; EElectrs+=energy;}
  else if(particle==positron) {nPositrs++; EPositrs+=energy;}
  else if(particle==proton)   {nProtons++; EProtons+=energy;}
  else if(particle==neutron)  {nNeutron++; ENeutron+=energy;}
  else if(particle==deuteron) {n_Deutrs++; E_Deutrs+=energy;}
  else if(particle==helium3)  {n_Heliu3++; E_Heliu3+=energy;}
  else if(particle==triton)   {n_Triton++; E_Triton+=energy;}
  else if(particle==alpha)    {n_Alphas++; E_Alphas+=energy;}
  else                        {n_Others++; E_Others+=energy;}
}

void Test23TrackingAction::PrintResult()
{
  if(n_Gammas)E_Gammas/=n_Gammas;
  if(n_NeutMu)E_NeutMu/=n_NeutMu;
  if(n_AnNuEl)E_AnNuEl/=n_AnNuEl;
  if(nElectrs)EElectrs/=nElectrs;
  if(nPositrs)EPositrs/=nPositrs;
  if(nProtons)EProtons/=nProtons;
  if(nNeutron)ENeutron/=nNeutron;
  if(n_Deutrs)E_Deutrs/=n_Deutrs;
  if(n_Heliu3)E_Heliu3/=n_Heliu3;
  if(n_Triton)E_Triton/=n_Triton;
  if(n_Alphas)E_Alphas/=n_Alphas;
  if(n_Others)E_Others/=n_Others;
  G4cout<<"Test23_Result: gam="<<n_Gammas<<"("<<E_Gammas<<"), e-="<<nElectrs<<"("<<EElectrs
        <<"), e+="<<nPositrs<<"("<<EPositrs<<"), p="<<nProtons<<"("<<EProtons<<"), n="
        <<nNeutron<<"("<<ENeutron<<"), d="<<n_Deutrs<<"("<<E_Deutrs<<"), t="<<n_Triton<<"("
        <<E_Triton<<"), He3="<<n_Heliu3<<"("<<E_Heliu3<<"), alph="<<n_Alphas<<"("<<E_Alphas
        <<"), nu="<<n_NeutMu<<"("<<E_NeutMu<<"), anu="<<n_AnNuEl<<"("<<E_AnNuEl
        <<"), others="<<n_Others<<"*"<<E_Others<<G4endl;
}

void Test23TrackingAction::ResetResult()
{
  n_Gammas=0; E_Gammas=0.;
  n_NeutMu=0; E_NeutMu=0.;
  n_AnNuEl=0; E_AnNuEl=0.;
  nElectrs=0; EElectrs=0.;
  nPositrs=0; EPositrs=0.;
  nProtons=0; EProtons=0.;
  nNeutron=0; ENeutron=0.;
  n_Deutrs=0; E_Deutrs=0.;
  n_Heliu3=0; E_Heliu3=0.;
  n_Triton=0; E_Triton=0.;
  n_Alphas=0; E_Alphas=0.;
  n_Others=0; E_Others=0.;
}
