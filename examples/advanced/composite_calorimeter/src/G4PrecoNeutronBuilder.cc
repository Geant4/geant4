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

G4PrecoNeutronBuilder::
G4PrecoNeutronBuilder() 
{
  theMin = 0;
  theMax = 150.*MeV;
  theModel = new G4PreCompoundModel(&theHandler);
}

G4PrecoNeutronBuilder::
~G4PrecoNeutronBuilder() {}

void G4PrecoNeutronBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4PrecoNeutronBuilder::
Build(G4HadronFissionProcess & aP)
{
}

void G4PrecoNeutronBuilder::
Build(G4HadronCaptureProcess & aP)
{
}

void G4PrecoNeutronBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theModel);
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(&theXSec);  
}

// 2002 by J.P. Wellisch
