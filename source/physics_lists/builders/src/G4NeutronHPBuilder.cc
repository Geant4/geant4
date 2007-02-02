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
#include "G4NeutronHPBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4NeutronHPBuilder::
G4NeutronHPBuilder() 
{
  theHPElastic = 0;
  theHPElasticData = 0;
  
  theHPFission = 0;
  theHPFissionData = 0;
  
  theHPCapture = 0;
  theHPCaptureData = 0;
  
  theHPInelastic = 0;
  theHPInelasticData = 0;

  theMin = 0;
  theIMin = theMin;
  theMax = 20*MeV;
  theIMax = theMax;

}

G4NeutronHPBuilder::
~G4NeutronHPBuilder() 
{
  delete theHPElasticData;
  delete theHPFissionData;
  delete theHPCaptureData;
  delete theHPInelasticData;
}

void G4NeutronHPBuilder::
Build(G4HadronElasticProcess * aP)
{
  if(theHPElastic==0) theHPElastic = new G4NeutronHPElastic;
  theHPElastic->SetMinEnergy(theMin);
  theHPElastic->SetMaxEnergy(theMax);
  if(theHPElasticData == 0) theHPElasticData = new G4NeutronHPElasticData;
  aP->AddDataSet(theHPElasticData);
  aP->RegisterMe(theHPElastic);
}

void G4NeutronHPBuilder::
Build(G4HadronFissionProcess * aP)
{
  if(theHPFission == 0) theHPFission = new G4NeutronHPFission;
  theHPFission->SetMinEnergy(theMin);
  theHPFission->SetMaxEnergy(theMax);
  if(theHPFissionData==0) theHPFissionData=new G4NeutronHPFissionData;
  aP->AddDataSet(theHPFissionData);
  aP->RegisterMe(theHPFission);
}

void G4NeutronHPBuilder::
Build(G4HadronCaptureProcess * aP)
{
  if(theHPCapture==0) theHPCapture = new G4NeutronHPCapture;
  theHPCapture->SetMinEnergy(theMin);
  theHPCapture->SetMaxEnergy(theMax);
  if(theHPCaptureData==0) theHPCaptureData = new G4NeutronHPCaptureData;
  aP->AddDataSet(theHPCaptureData);
  aP->RegisterMe(theHPCapture);
}

void G4NeutronHPBuilder::
Build(G4NeutronInelasticProcess * aP)
{
  if(theHPInelastic==0) theHPInelastic = new G4NeutronHPInelastic;
  theHPInelastic->SetMinEnergy(theIMin);
  theHPInelastic->SetMaxEnergy(theIMax);
  if(theHPInelasticData==0) theHPInelasticData = new G4NeutronHPInelasticData;
  aP->AddDataSet(theHPInelasticData);
  aP->RegisterMe(theHPInelastic);
}
// 2002 by J.P. Wellisch
