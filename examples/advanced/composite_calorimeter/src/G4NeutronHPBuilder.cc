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

G4NeutronHPBuilder::
G4NeutronHPBuilder() 
{
  theHPElastic = NULL;
  theHPElasticData = NULL;
  
  theHPFission = NULL;
  theHPFissionData = NULL;
  
  theHPCapture = NULL;
  theHPCaptureData = NULL;
  
  theHPInelastic = NULL;
  theHPInelasticData = NULL;
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
Build(G4HadronElasticProcess & aP)
{
  if(theHPElastic==NULL) theHPElastic = new G4NeutronHPElastic;
  if(theHPElasticData == NULL) theHPElasticData = new G4NeutronHPElasticData;
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(theHPElasticData);
  aP.RegisterMe(theHPElastic);
}

void G4NeutronHPBuilder::
Build(G4HadronFissionProcess & aP)
{
  if(theHPFission == NULL) theHPFission = new G4NeutronHPFission;
  if(theHPFissionData==NULL) theHPFissionData=new G4NeutronHPFissionData;
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(theHPFissionData);
  aP.RegisterMe(theHPFission);
}

void G4NeutronHPBuilder::
Build(G4HadronCaptureProcess & aP)
{
  if(theHPCapture==NULL) theHPCapture = new G4NeutronHPCapture;
  if(theHPCaptureData==NULL) theHPCaptureData = new G4NeutronHPCaptureData;
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(theHPCaptureData);
  aP.RegisterMe(theHPCapture);
}

void G4NeutronHPBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  if(theHPInelastic==NULL) theHPInelastic = new G4NeutronHPInelastic;
  if(theHPInelasticData==NULL) theHPInelasticData = new G4NeutronHPInelasticData;
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(theHPInelasticData);
  aP.RegisterMe(theHPInelastic);
}
// 2002 by J.P. Wellisch
