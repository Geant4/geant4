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
// $Id: G4PiKBuilder.cc,v 1.2 2005-11-25 15:38:50 gunter Exp $
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
G4PiKBuilder(): wasActivated(false) {}

G4PiKBuilder::
~G4PiKBuilder(){
  if(wasActivated)
  {
    G4ProcessManager * theProcMan;
    theProcMan = G4PionPlus::PionPlus()->GetProcessManager();
    if(theProcMan) theProcMan->RemoveProcess(thePionPlusElasticProcess);
    if(theProcMan) theProcMan->RemoveProcess(thePionPlusInelastic);
    theProcMan = G4PionMinus::PionMinus()->GetProcessManager();
    if(theProcMan) theProcMan->RemoveProcess(thePionMinusElasticProcess);
    if(theProcMan) theProcMan->RemoveProcess(thePionMinusInelastic);
    theProcMan = G4KaonPlus::KaonPlus()->GetProcessManager();
    if(theProcMan) theProcMan->RemoveProcess(theKaonPlusElasticProcess);
    if(theProcMan) theProcMan->RemoveProcess(theKaonPlusInelastic);
    theProcMan = G4KaonMinus::KaonMinus()->GetProcessManager();
    if(theProcMan) theProcMan->RemoveProcess(theKaonMinusElasticProcess);
    if(theProcMan) theProcMan->RemoveProcess(theKaonMinusInelastic);
    theProcMan = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
    if(theProcMan) theProcMan->RemoveProcess(theKaonZeroLElasticProcess);
    if(theProcMan) theProcMan->RemoveProcess(theKaonZeroLInelastic);
    theProcMan = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();
    if(theProcMan) theProcMan->RemoveProcess(theKaonZeroSElasticProcess);
    if(theProcMan) theProcMan->RemoveProcess(theKaonZeroSInelastic);
  }
}

void G4PiKBuilder::
Build()
{
  wasActivated = true;
  std::vector<G4VPiKBuilder *>::iterator i;
  for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
  {
    (*i)->Build(thePionPlusElasticProcess=new G4HadronElasticProcess);
    (*i)->Build(thePionMinusElasticProcess=new G4HadronElasticProcess);
    (*i)->Build(theKaonPlusElasticProcess=new G4HadronElasticProcess);
    (*i)->Build(theKaonMinusElasticProcess=new G4HadronElasticProcess);
    (*i)->Build(theKaonZeroLElasticProcess=new G4HadronElasticProcess);
    (*i)->Build(theKaonZeroSElasticProcess=new G4HadronElasticProcess);

    (*i)->Build(thePionPlusInelastic=new G4PionPlusInelasticProcess);
    (*i)->Build(thePionMinusInelastic=new G4PionMinusInelasticProcess);
    (*i)->Build(theKaonPlusInelastic=new G4KaonPlusInelasticProcess);
    (*i)->Build(theKaonMinusInelastic=new G4KaonMinusInelasticProcess);
    (*i)->Build(theKaonZeroLInelastic=new G4KaonZeroLInelasticProcess);
    (*i)->Build(theKaonZeroSInelastic=new G4KaonZeroSInelasticProcess);
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
