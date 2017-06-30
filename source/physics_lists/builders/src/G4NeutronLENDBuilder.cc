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
#include "G4NeutronLENDBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4NeutronLENDBuilder::
G4NeutronLENDBuilder( G4String eva ) 
{
  theLENDElastic = 0;
  theLENDElasticCrossSection = 0;
  
  theLENDFission = 0;
  theLENDFissionCrossSection = 0;
  
  theLENDCapture = 0;
  theLENDCaptureCrossSection = 0;
  
  theLENDInelastic = 0;
  theLENDInelasticCrossSection = 0;

  theMin = 0;
  theIMin = theMin;
  theMax = 20*MeV;
  theIMax = theMax;
  evaluation = eva;

}
void G4NeutronLENDBuilder::
Build(G4HadronElasticProcess * aP)
{
  if(theLENDElastic==0) theLENDElastic = new G4LENDElastic( G4Neutron::Neutron() );
  theLENDElastic->SetMinEnergy(theMin);
  theLENDElastic->SetMaxEnergy(theMax);

  if ( evaluation != "" ) theLENDElastic->ChangeDefaultEvaluation( evaluation );
  //theLENDElastic->AllowNaturalAbundanceTarget();
  theLENDElastic->AllowAnyCandidateTarget();
  if(theLENDElasticCrossSection == 0) theLENDElasticCrossSection = new G4LENDElasticCrossSection( G4Neutron::Neutron() );
  if ( evaluation != "" ) theLENDElasticCrossSection->ChangeDefaultEvaluation( evaluation );
  //theLENDElasticCrossSection->AllowNaturalAbundanceTarget();
  theLENDElasticCrossSection->AllowAnyCandidateTarget();
  aP->AddDataSet(theLENDElasticCrossSection);
  aP->RegisterMe(theLENDElastic);
}

void G4NeutronLENDBuilder::
Build(G4HadronFissionProcess * aP)
{
  if(theLENDFission == 0) theLENDFission = new G4LENDFission( G4Neutron::Neutron() );
  theLENDFission->SetMinEnergy(theMin);
  theLENDFission->SetMaxEnergy(theMax);
  if ( evaluation != "" ) theLENDFission->ChangeDefaultEvaluation( evaluation );
  //theLENDFission->AllowNaturalAbundanceTarget();
  theLENDFission->AllowAnyCandidateTarget();
  if(theLENDFissionCrossSection==0) theLENDFissionCrossSection=new G4LENDFissionCrossSection( G4Neutron::Neutron() );
  if ( evaluation != "" ) theLENDFissionCrossSection->ChangeDefaultEvaluation( evaluation );
  //theLENDFissionCrossSection->AllowNaturalAbundanceTarget();
  theLENDFissionCrossSection->AllowAnyCandidateTarget();
  aP->AddDataSet(theLENDFissionCrossSection);
  aP->RegisterMe(theLENDFission);
}

void G4NeutronLENDBuilder::
Build(G4HadronCaptureProcess * aP)
{
  if(theLENDCapture==0) theLENDCapture = new G4LENDCapture( G4Neutron::Neutron() );
  theLENDCapture->SetMinEnergy(theMin);
  theLENDCapture->SetMaxEnergy(theMax);
  if ( evaluation != "" ) theLENDCapture->ChangeDefaultEvaluation( evaluation );
  //theLENDCapture->AllowNaturalAbundanceTarget();
  theLENDCapture->AllowAnyCandidateTarget();
  if(theLENDCaptureCrossSection==0) theLENDCaptureCrossSection = new G4LENDCaptureCrossSection( G4Neutron::Neutron() );
  if ( evaluation != "" ) theLENDCaptureCrossSection->ChangeDefaultEvaluation( evaluation );
  //theLENDCaptureCrossSection->AllowNaturalAbundanceTarget();
  theLENDCaptureCrossSection->AllowAnyCandidateTarget();
  aP->AddDataSet(theLENDCaptureCrossSection);
  aP->RegisterMe(theLENDCapture);
}

void G4NeutronLENDBuilder::
Build(G4NeutronInelasticProcess * aP)
{
  if(theLENDInelastic==0) theLENDInelastic = new G4LENDInelastic( G4Neutron::Neutron() );
  theLENDInelastic->SetMinEnergy(theIMin);
  theLENDInelastic->SetMaxEnergy(theIMax);
  if ( evaluation != "" ) theLENDInelastic->ChangeDefaultEvaluation( evaluation );
  //theLENDInelastic->AllowNaturalAbundanceTarget();
  theLENDInelastic->AllowAnyCandidateTarget();
  if(theLENDInelasticCrossSection==0) theLENDInelasticCrossSection = new G4LENDInelasticCrossSection( G4Neutron::Neutron() );
  if ( evaluation != "" ) theLENDInelasticCrossSection->ChangeDefaultEvaluation( evaluation );
  //theLENDInelasticCrossSection->AllowNaturalAbundanceTarget();
  theLENDInelasticCrossSection->AllowAnyCandidateTarget();
  aP->AddDataSet(theLENDInelasticCrossSection);
  aP->RegisterMe(theLENDInelastic);
}
