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
#include "G4INCLXXProtonBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4PreCompoundModel.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronicInteractionRegistry.hh"

G4INCLXXProtonBuilder::
G4INCLXXProtonBuilder() 
{
  thePreCompoundMin = 0;
  thePreCompoundMax = 2*MeV;
  theMin = 1.0*MeV;
  theMax = 3.0*GeV;
  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  thePreCompoundModel = static_cast<G4VPreCompoundModel*>(p);
  if(!thePreCompoundModel) { thePreCompoundModel = new G4PreCompoundModel(); }
  theModel = new G4INCLXXInterface(thePreCompoundModel);
}

void G4INCLXXProtonBuilder::
Build(G4ProtonInelasticProcess * aP)
{
  thePreCompoundModel->SetMinEnergy(thePreCompoundMin);
  thePreCompoundModel->SetMaxEnergy(thePreCompoundMax);
  aP->RegisterMe(thePreCompoundModel);
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theModel);
  aP->AddDataSet(new G4BGGNucleonInelasticXS(G4Proton::Proton()));
}

// 2011 by P. Kaitaniemi

