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

G4PiKBuilder::
G4PiKBuilder() {}

G4PiKBuilder::
~G4PiKBuilder() {}

void G4PiKBuilder::
Build()
{
  G4std::vector<G4VPiKBuilder *>::iterator i;
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
  theProcMan->AddDiscreteProcess(&thePionPlusElasticProcess);
  theProcMan->AddDiscreteProcess(&thePionPlusInelastic);
  theProcMan = G4PionMinus::PionMinus()->GetProcessManager();
  theProcMan->AddDiscreteProcess(&thePionMinusElasticProcess);
  theProcMan->AddDiscreteProcess(&thePionMinusInelastic);
  theProcMan = G4KaonPlus::KaonPlus()->GetProcessManager();
  theProcMan->AddDiscreteProcess(&theKaonPlusElasticProcess);
  theProcMan->AddDiscreteProcess(&theKaonPlusInelastic);
  theProcMan = G4KaonMinus::KaonMinus()->GetProcessManager();
  theProcMan->AddDiscreteProcess(&theKaonMinusElasticProcess);
  theProcMan->AddDiscreteProcess(&theKaonMinusInelastic);
  theProcMan = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
  theProcMan->AddDiscreteProcess(&theKaonZeroLElasticProcess);
  theProcMan->AddDiscreteProcess(&theKaonZeroLInelastic);
  theProcMan = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();
  theProcMan->AddDiscreteProcess(&theKaonZeroSElasticProcess);
  theProcMan->AddDiscreteProcess(&theKaonZeroSInelastic);
}
// 2002 by J.P. Wellisch
