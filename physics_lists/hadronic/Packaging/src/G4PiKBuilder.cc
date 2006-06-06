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
// $Id: G4PiKBuilder.cc,v 1.5 2006-06-06 16:47:45 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4PiKBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 16.11.2005 G.Folger: don't  keep processes as data members, but new these
//
//----------------------------------------------------------------------------
//
#include "G4PiKBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4PiKBuilder::
G4PiKBuilder(): wasActivated(false) 
{  
  thePionPlusElasticProcess=new G4HadronElasticProcess;
  thePionMinusElasticProcess=new G4HadronElasticProcess;
  theKaonPlusElasticProcess=new G4HadronElasticProcess;
  theKaonMinusElasticProcess=new G4HadronElasticProcess;
  theKaonZeroLElasticProcess=new G4HadronElasticProcess;
  theKaonZeroSElasticProcess=new G4HadronElasticProcess;

  thePionPlusInelastic=new G4PionPlusInelasticProcess;
  thePionMinusInelastic=new G4PionMinusInelasticProcess;
  theKaonPlusInelastic=new G4KaonPlusInelasticProcess;
  theKaonMinusInelastic=new G4KaonMinusInelasticProcess;
  theKaonZeroLInelastic=new G4KaonZeroLInelasticProcess;
  theKaonZeroSInelastic=new G4KaonZeroSInelasticProcess;
}

G4PiKBuilder::
~G4PiKBuilder(){
  delete thePionPlusElasticProcess;
  delete thePionMinusElasticProcess;
  delete theKaonPlusElasticProcess;
  delete theKaonMinusElasticProcess;
  delete theKaonZeroLElasticProcess;
  delete theKaonZeroSElasticProcess;

  delete thePionPlusInelastic;
  delete thePionMinusInelastic;
  delete theKaonPlusInelastic;
  delete theKaonMinusInelastic;
  delete theKaonZeroLInelastic;
  delete theKaonZeroSInelastic;
}

void G4PiKBuilder::
Build()
{
  wasActivated = true;

  std::vector<G4VPiKBuilder *>::iterator i;
  for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
  {
    (*i)->Build(thePionPlusElasticProcess);
    (*i)->Build(thePionMinusElasticProcess);
    (*i)->Build(theKaonPlusElasticProcess);
    (*i)->Build(theKaonMinusElasticProcess);
    (*i)->Build(theKaonZeroLElasticProcess);
    (*i)->Build(theKaonZeroSElasticProcess);

    (*i)->Build(thePionPlusInelastic);
    (*i)->Build(thePionMinusInelastic);
    (*i)->Build(theKaonPlusInelastic);
    (*i)->Build(theKaonMinusInelastic);
    (*i)->Build(theKaonZeroLInelastic);
    (*i)->Build(theKaonZeroSInelastic);
  }
  G4ProcessManager * theProcMan;
  theProcMan = G4PionPlus::PionPlus()->GetProcessManager();
  theProcMan->AddDiscreteProcess(thePionPlusElasticProcess);
  theProcMan->AddDiscreteProcess(thePionPlusInelastic);
  
  theProcMan = G4PionMinus::PionMinus()->GetProcessManager();
  theProcMan->AddDiscreteProcess(thePionMinusElasticProcess);
  theProcMan->AddDiscreteProcess(thePionMinusInelastic);
  
  theProcMan = G4KaonPlus::KaonPlus()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theKaonPlusElasticProcess);
  theProcMan->AddDiscreteProcess(theKaonPlusInelastic);
  
  theProcMan = G4KaonMinus::KaonMinus()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theKaonMinusElasticProcess);
  theProcMan->AddDiscreteProcess(theKaonMinusInelastic);
  
  theProcMan = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theKaonZeroLElasticProcess);
  theProcMan->AddDiscreteProcess(theKaonZeroLInelastic);
  
  theProcMan = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();
  theProcMan->AddDiscreteProcess(theKaonZeroSElasticProcess);
  theProcMan->AddDiscreteProcess(theKaonZeroSInelastic);
}
// 2002 by J.P. Wellisch
