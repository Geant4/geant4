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
#include "G4BinaryProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4BinaryProtonBuilder::
G4BinaryProtonBuilder() 
{
  theModel = new G4BinaryCascade;
  theMin = 0;
  theMax = 9.9*GeV;
}

void G4BinaryProtonBuilder::
Build(G4ProtonInelasticProcess & aP)
{
  aP.AddDataSet(&theXSec);  
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theModel);
}

G4BinaryProtonBuilder::
~G4BinaryProtonBuilder() {}

void G4BinaryProtonBuilder::
Build(G4HadronElasticProcess & )
{
}

// 2002 by J.P. Wellisch
