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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4KaonBuilder
//
// Author: 2010 G.Folger
//  devired from G4PiKBuilder
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "G4KaonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4KaonBuilder::
G4KaonBuilder(): wasActivated(false) 
{  
  theKaonPlusInelastic=new G4KaonPlusInelasticProcess;
  theKaonMinusInelastic=new G4KaonMinusInelasticProcess;
  theKaonZeroLInelastic=new G4KaonZeroLInelasticProcess;
  theKaonZeroSInelastic=new G4KaonZeroSInelasticProcess;
}

G4KaonBuilder::
~G4KaonBuilder(){
  delete theKaonPlusInelastic;
  delete theKaonMinusInelastic;
  delete theKaonZeroLInelastic;
  delete theKaonZeroSInelastic;
}

void G4KaonBuilder::
Build()
{
  wasActivated = true;

  std::vector<G4VKaonBuilder *>::iterator i;
  for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
  {
    (*i)->Build(theKaonPlusInelastic);
    (*i)->Build(theKaonMinusInelastic);
    (*i)->Build(theKaonZeroLInelastic);
    (*i)->Build(theKaonZeroSInelastic);
  }
  G4ProcessManager * theProcMan;
 
  theProcMan = G4KaonPlus::KaonPlus()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theKaonPlusInelastic);
  
  theProcMan = G4KaonMinus::KaonMinus()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theKaonMinusInelastic);
  
  theProcMan = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theKaonZeroLInelastic);
  
  theProcMan = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theKaonZeroSInelastic);
}
