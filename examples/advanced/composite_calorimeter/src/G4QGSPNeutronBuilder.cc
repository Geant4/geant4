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
#include "G4ProcessManager.hh"

G4QGSPNeutronBuilder::
G4QGSPNeutronBuilder() 
{
  theMin = 15*GeV;
  theModel = new G4TheoFSGenerator;
  theCascade = new G4GeneratorPrecompoundInterface;
  thePreEquilib = new G4PreCompoundModel(&theHandler);
  theCascade->SetDeExcitation(thePreEquilib);  
  theModel->SetTransport(theCascade);
  theModel->SetHighEnergyGenerator(&theStringModel);
  theStringDecay = new G4ExcitedStringDecay(&theFragmentation);
  theStringModel.SetFragmentationModel(theStringDecay);
}

G4QGSPNeutronBuilder::
~G4QGSPNeutronBuilder() 
{
  delete theStringDecay;
}

void G4QGSPNeutronBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4QGSPNeutronBuilder::
Build(G4HadronFissionProcess & aP)
{
}

void G4QGSPNeutronBuilder::
Build(G4HadronCaptureProcess & aP)
{
}

void G4QGSPNeutronBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
  aP.RegisterMe(theModel);
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(&theXSec);  
}

// 2002 by J.P. Wellisch
