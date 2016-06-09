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
  if(theHPElasticData == 0) theHPElasticData = new G4NeutronHPElasticData;
  aP->AddDataSet(theHPElasticData);
  aP->RegisterMe(theHPElastic);
}

void G4NeutronHPBuilder::
Build(G4HadronFissionProcess * aP)
{
  if(theHPFission == 0) theHPFission = new G4NeutronHPFission;
  if(theHPFissionData==0) theHPFissionData=new G4NeutronHPFissionData;
  aP->AddDataSet(theHPFissionData);
  aP->RegisterMe(theHPFission);
}

void G4NeutronHPBuilder::
Build(G4HadronCaptureProcess * aP)
{
  if(theHPCapture==0) theHPCapture = new G4NeutronHPCapture;
  if(theHPCaptureData==0) theHPCaptureData = new G4NeutronHPCaptureData;
  aP->AddDataSet(theHPCaptureData);
  aP->RegisterMe(theHPCapture);
}

void G4NeutronHPBuilder::
Build(G4NeutronInelasticProcess * aP)
{
  if(theHPInelastic==0) theHPInelastic = new G4NeutronHPInelastic;
  if(theHPInelasticData==0) theHPInelasticData = new G4NeutronHPInelasticData;
  aP->AddDataSet(theHPInelasticData);
  aP->RegisterMe(theHPInelastic);
}
// 2002 by J.P. Wellisch
