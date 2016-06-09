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
#include "G4BertiniNeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4BertiniNeutronBuilder::
G4BertiniNeutronBuilder() 
{
  theMin = 0;
  theMax = 5*GeV;
  theModel = new G4CascadeInterface;
}

G4BertiniNeutronBuilder::
~G4BertiniNeutronBuilder() {}

void G4BertiniNeutronBuilder::
Build(G4HadronElasticProcess & )
{
}

void G4BertiniNeutronBuilder::
Build(G4HadronFissionProcess & )
{
}

void G4BertiniNeutronBuilder::
Build(G4HadronCaptureProcess & )
{
}

void G4BertiniNeutronBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theModel);
  aP.AddDataSet(&theXSec);  
}

// 2002 by J.P. Wellisch
