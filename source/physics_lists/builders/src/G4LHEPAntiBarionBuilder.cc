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
// ClassName:   G4LHEPAntiBarionBuilder
//
// Author: 2011 J. Apostolakis
//
// Modified:
// 16.11.2005 G.Folger: don't  keep processes as data members, but new these
// 13.06.2006 G.Folger: (re)move elastic scatterring 
//
//----------------------------------------------------------------------------
//
#include "G4LHEPAntiBarionBuilder.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"


G4LHEPAntiBarionBuilder::G4LHEPAntiBarionBuilder():  
 theAntiProtonInelastic(0), theLEAntiProtonModel(0), theHEAntiProtonModel(0),
 theAntiNeutronInelastic(0), theLEAntiNeutronModel(0), theHEAntiNeutronModel(0),
 wasActivated(false)
{}

G4LHEPAntiBarionBuilder::~G4LHEPAntiBarionBuilder()
{}

void G4LHEPAntiBarionBuilder::Build()
{
  G4ProcessManager * aProcMan = 0;
  wasActivated = true;
  
  // anti-Proton
  theAntiProtonInelastic = new G4AntiProtonInelasticProcess();
  aProcMan = G4AntiProton::AntiProton()->GetProcessManager();
  theLEAntiProtonModel = new G4LEAntiProtonInelastic();
  theHEAntiProtonModel = new G4HEAntiProtonInelastic();
  theHEAntiProtonModel->SetMaxEnergy(100*TeV);
  theAntiProtonInelastic->RegisterMe(theLEAntiProtonModel);
  theAntiProtonInelastic->RegisterMe(theHEAntiProtonModel);
  aProcMan->AddDiscreteProcess(theAntiProtonInelastic);

  // AntiNeutron
  theAntiNeutronInelastic = new G4AntiNeutronInelasticProcess();
  aProcMan = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  theLEAntiNeutronModel = new G4LEAntiNeutronInelastic();
  theHEAntiNeutronModel = new G4HEAntiNeutronInelastic();
  theHEAntiNeutronModel->SetMaxEnergy(100*TeV);
  theAntiNeutronInelastic->RegisterMe(theLEAntiNeutronModel);
  theAntiNeutronInelastic->RegisterMe(theHEAntiNeutronModel);
  aProcMan->AddDiscreteProcess(theAntiNeutronInelastic);
}
