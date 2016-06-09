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
#include "G4LEADPiKBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4LEADPiKBuilder::
G4LEADPiKBuilder() 
{
  theMin = 0;
  theModel = new G4Mars5GeV;
}

G4LEADPiKBuilder::
~G4LEADPiKBuilder() {}

void G4LEADPiKBuilder::
Build(G4HadronElasticProcess & ) {}

void G4LEADPiKBuilder::
Build(G4PionPlusInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  aP.AddDataSet(&thePiData);
  theModel->SetMinEnergy(theMin);
}

void G4LEADPiKBuilder::
Build(G4PionMinusInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  aP.AddDataSet(&thePiData);
  theModel->SetMinEnergy(theMin);
}

void G4LEADPiKBuilder::
Build(G4KaonPlusInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  theModel->SetMinEnergy(theMin);
}

void G4LEADPiKBuilder::
Build(G4KaonMinusInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  theModel->SetMinEnergy(theMin);
}

void G4LEADPiKBuilder::
Build(G4KaonZeroLInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  theModel->SetMinEnergy(theMin);
}

void G4LEADPiKBuilder::
Build(G4KaonZeroSInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  theModel->SetMinEnergy(theMin);
}

// 2002 by J.P. Wellisch
