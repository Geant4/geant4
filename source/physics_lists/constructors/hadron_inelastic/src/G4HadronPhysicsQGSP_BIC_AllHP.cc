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
// $Id: G4HadronPhysicsQGSP_BIC_AllHP.cc 73040 2013-08-15 09:36:57Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsQGSP_BIC_AllHP
//
// Author: 2013 P. Arce
//
// based on G4HadronPhysicsQGSP_BIC_HP
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4SystemOfUnits.hh"
#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"
#include "G4ProtonPHPBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "G4QGSPProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"
#include "G4ProtonPHPBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
#include "G4BinaryNeutronBuilder.hh"
#include "G4NeutronPHPBuilder.hh"




#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_BIC_AllHP);

G4HadronPhysicsQGSP_BIC_AllHP::G4HadronPhysicsQGSP_BIC_AllHP(G4int)
    :  G4HadronPhysicsQGSP_BIC_AllHP("hInelastic QGSP_BIC_HP")
{}

G4HadronPhysicsQGSP_BIC_AllHP::G4HadronPhysicsQGSP_BIC_AllHP(const G4String& name, G4bool /* quasiElastic */)
    :  G4HadronPhysicsQGSP_BIC(name)  
{
   minBIC_neutron = 19.9*MeV;
   maxHP_neutron = 20.*MeV;
   minBIC_proton = 200.*MeV;
   maxHP_proton = 200.*MeV; 
}

void G4HadronPhysicsQGSP_BIC_AllHP::Neutron()
{
  auto neu = new G4NeutronBuilder( true ); // Fission on
  AddBuilder(neu);
  auto qgs = new G4QGSPNeutronBuilder(QuasiElasticQGS);
  AddBuilder(qgs);
  qgs->SetMinEnergy(minQGSP_neutron);
  neu->RegisterMe(qgs);
  auto ftf = new G4FTFPNeutronBuilder(QuasiElasticFTF);
  AddBuilder(ftf);
  ftf->SetMinEnergy(minFTFP_neutron);
  ftf->SetMaxEnergy(maxFTFP_neutron);
  neu->RegisterMe(ftf);
  auto bic = new G4BinaryNeutronBuilder;
  AddBuilder(bic);
  bic->SetMinEnergy(minBIC_neutron);
  bic->SetMaxEnergy(maxBIC_neutron);
  neu->RegisterMe(bic);
  auto hp = new G4NeutronPHPBuilder;
  AddBuilder(hp);
  hp->SetMaxEnergy(maxHP_neutron);
  neu->RegisterMe(hp);
  neu->Build();
}

void G4HadronPhysicsQGSP_BIC_AllHP::Proton()
{
  auto pro = new G4ProtonBuilder;
  AddBuilder(pro);
  auto qgs = new G4QGSPProtonBuilder(QuasiElasticQGS);
  AddBuilder(qgs);
  qgs->SetMinEnergy(minQGSP_proton);
  pro->RegisterMe(qgs);
  auto ftf = new G4FTFPProtonBuilder(QuasiElasticFTF);
  AddBuilder(ftf);
  ftf->SetMinEnergy(minFTFP_proton);
  ftf->SetMaxEnergy(maxFTFP_proton);
  pro->RegisterMe(ftf);
  auto bic = new G4BinaryProtonBuilder;
  AddBuilder(bic);
  bic->SetMaxEnergy(maxBIC_proton);
  bic->SetMinEnergy(minBIC_proton);
  pro->RegisterMe(bic);
  auto hp = new G4ProtonPHPBuilder;
  AddBuilder(hp);
  hp->SetMaxEnergy(maxHP_proton);
  pro->RegisterMe(hp);
  pro->Build();
}

