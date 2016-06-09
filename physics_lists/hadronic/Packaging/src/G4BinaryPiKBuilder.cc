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
#include "G4BinaryPiKBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4BinaryPiKBuilder::
G4BinaryPiKBuilder() 
{
  theMin = 0*GeV;
  theMax = 1.3*GeV;
  theModel = new G4BinaryCascade;
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax); 
}

G4BinaryPiKBuilder::
~G4BinaryPiKBuilder() {}

void G4BinaryPiKBuilder::
Build(G4HadronElasticProcess * ) {}

void G4BinaryPiKBuilder::
Build(G4PionPlusInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theModel);
  aP->AddDataSet(&thePiData);
}

void G4BinaryPiKBuilder::
Build(G4PionMinusInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theModel);
  aP->AddDataSet(&thePiData);
}

void G4BinaryPiKBuilder::
Build(G4KaonPlusInelasticProcess * )
{
}

void G4BinaryPiKBuilder::
Build(G4KaonMinusInelasticProcess * )
{
}

void G4BinaryPiKBuilder::
Build(G4KaonZeroLInelasticProcess * )
{
}

void G4BinaryPiKBuilder::
Build(G4KaonZeroSInelasticProcess * )
{
}

// 2002 by J.P. Wellisch
