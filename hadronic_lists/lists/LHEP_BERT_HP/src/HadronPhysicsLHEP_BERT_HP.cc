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
#include "HadronPhysicsLHEP_BERT_HP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsLHEP_BERT_HP::HadronPhysicsLHEP_BERT_HP(const G4String& name)
                    :  G4VPhysicsConstructor(name) 
{
  theNeutrons.RegisterMe(&theLHEPNeutron);
  theNeutrons.RegisterMe(&theBertiniNeutron);
  theNeutrons.RegisterMe(&theHPNeutron);
  theLHEPNeutron.SetMinEnergy(19.9*MeV);
  theLHEPNeutron.SetMinInelasticEnergy(2.8*GeV);
  theBertiniNeutron.SetMaxEnergy(3.2*GeV);
  theBertiniNeutron.SetMinEnergy(19.9*MeV);

  thePro.RegisterMe(&theLHEPPro);
  thePro.RegisterMe(&theBertiniPro);
  theLHEPPro.SetMinEnergy(2.8*GeV);
  theBertiniPro.SetMaxEnergy(3.2*GeV);
  
  thePiK.RegisterMe(&theLHEPPiK);
  
}

HadronPhysicsLHEP_BERT_HP::~HadronPhysicsLHEP_BERT_HP() {}

void HadronPhysicsLHEP_BERT_HP::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsLHEP_BERT_HP::ConstructProcess()
{
  theNeutrons.Build();
  thePro.Build();
  thePiK.Build();
  theMiscLHEP.Build();
  theStoppingHadron.Build();
  theHadronQED.Build();
}
// 2002 by J.P. Wellisch
