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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4LEPNeutronBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 16.11.2005 G.Folger: don't  keep processes as data members, but new these
// 13.06.2006 G.Folger: (re)move elastic scatterring 
//
//----------------------------------------------------------------------------
//
#include "G4LEPNeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4LEPNeutronBuilder::G4LEPNeutronBuilder() : theLENeutronModel(0),
					     theNeutronFissionModel(0),
					     theNeutronCaptureModel(0) {
  theMin = 0;
  theIMin = theMin;
  theMax = 20*TeV;
  theIMax = 55*GeV;
}

G4LEPNeutronBuilder::~G4LEPNeutronBuilder() 
{
  if (theNeutronFissionModel) delete theNeutronFissionModel;
  if (theNeutronCaptureModel) delete theNeutronCaptureModel;
  if (theLENeutronModel)      delete theLENeutronModel;
}

void G4LEPNeutronBuilder::
Build(G4HadronElasticProcess *)
{
     G4cout << "Info - G4LEPNeutronBuilder::Build() not adding elastic" << G4endl;
}

void G4LEPNeutronBuilder::
Build(G4HadronFissionProcess * aP)
{
  theNeutronFissionModel = new G4LFission();
  theNeutronFissionModel->SetMinEnergy(theMin);
  theNeutronFissionModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theNeutronFissionModel);
}

void G4LEPNeutronBuilder::
Build(G4HadronCaptureProcess * aP)
{
  theNeutronCaptureModel = new G4LCapture();
  theNeutronCaptureModel->SetMinEnergy(theMin);
  theNeutronCaptureModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theNeutronCaptureModel);
}

void G4LEPNeutronBuilder::
Build(G4NeutronInelasticProcess * aP)
{
  if ( theIMax > 1.*eV ) 
  {
     theLENeutronModel = new G4LENeutronInelastic();
     theLENeutronModel->SetMinEnergy(theIMin);
     theLENeutronModel->SetMaxEnergy(theIMax);
     aP->RegisterMe(theLENeutronModel);
  }
}

// 2002 by J.P. Wellisch
