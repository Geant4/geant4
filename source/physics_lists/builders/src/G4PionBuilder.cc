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
//
//---------------------------------------------------------------------------
//
// ClassName:   G4PionBuilder
//
// Author: 2010 G.Folger
//  devired from G4PiKBuilder
//
// Modified:
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#include "G4PionBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4PionBuilder::
G4PionBuilder(): wasActivated(false) 
{  
  thePionPlusInelastic=new G4PionPlusInelasticProcess;
  thePionMinusInelastic=new G4PionMinusInelasticProcess;
}

void G4PionBuilder::
Build()
{
  wasActivated = true;

  std::vector<G4VPionBuilder *>::iterator i;
  for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
  {
    (*i)->Build(thePionPlusInelastic);
    (*i)->Build(thePionMinusInelastic);
  }
  G4ProcessManager * theProcMan;

  theProcMan = G4PionPlus::PionPlus()->GetProcessManager();
  theProcMan->AddDiscreteProcess(thePionPlusInelastic);
  
  theProcMan = G4PionMinus::PionMinus()->GetProcessManager();
  theProcMan->AddDiscreteProcess(thePionMinusInelastic);
}

void G4PionBuilder::RegisterMe(G4PhysicsBuilderInterface* aB) {
  auto bld = dynamic_cast<G4VPionBuilder*>(aB);
  if ( bld != nullptr ) {
      theModelCollections.push_back(bld);
  } else {
      G4PhysicsBuilderInterface::RegisterMe(aB);
  }
}

