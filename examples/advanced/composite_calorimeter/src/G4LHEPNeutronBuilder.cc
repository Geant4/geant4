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
#include "G4LHEPNeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4LHEPNeutronBuilder::
G4LHEPNeutronBuilder() 
{
  theMin = 0;
  theIMin = 0;
}

G4LHEPNeutronBuilder::
~G4LHEPNeutronBuilder() {}

void G4LHEPNeutronBuilder::
Build(G4HadronElasticProcess & aP)
{
  theElasticModel = new G4LElastic();
  theElasticModel->SetMinEnergy(theMin);
  aP.RegisterMe(theElasticModel);
}

void G4LHEPNeutronBuilder::
Build(G4HadronFissionProcess & aP)
{
  theNeutronFissionModel = new G4LFission();
  theNeutronFissionModel->SetMinEnergy(theMin);
  aP.RegisterMe(theNeutronFissionModel);
}

void G4LHEPNeutronBuilder::
Build(G4HadronCaptureProcess & aP)
{
  theNeutronCaptureModel = new G4LCapture();
  theNeutronCaptureModel->SetMinEnergy(theMin);
  aP.RegisterMe(theNeutronCaptureModel);
}

void G4LHEPNeutronBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  theLENeutronModel = new G4LENeutronInelastic();
  theHENeutronModel = new G4HENeutronInelastic();
  theLENeutronModel->SetMinEnergy(theIMin);
  theLENeutronModel->SetMaxEnergy(55*GeV);
  theHENeutronModel->SetMinEnergy(25*GeV);
  aP.RegisterMe(theLENeutronModel);
  aP.RegisterMe(theHENeutronModel);
}

// 2002 by J.P. Wellisch
