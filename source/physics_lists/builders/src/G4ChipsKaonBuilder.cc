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
// Author: 2010 G.Folger
//  devired from G4LEPPiKBuilder
//
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "G4ChipsKaonBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4MesonConstructor.hh"
#include "G4ProcessManager.hh"

G4ChipsKaonBuilder::G4ChipsKaonBuilder(G4int ver)
  : verb(ver)
{
  theInelastic = new G4QInelastic();
}

G4ChipsKaonBuilder::
~G4ChipsKaonBuilder() 
{
  delete theInelastic;
}

void G4ChipsKaonBuilder::Build()
{
     static G4bool onceOnly(true);
     if ( onceOnly )
     {
	if (verb > 0 ) 
	  {G4cout << "Info - G4ChipsKaonBuilder::Build() not adding elastic" <<
                      G4endl;}

	attachProcess(G4KaonPlus::KaonPlus());
	attachProcess(G4KaonMinus::KaonMinus());
	attachProcess(G4KaonZeroShort::KaonZeroShort());
	attachProcess(G4KaonZeroLong::KaonZeroLong());
        onceOnly=false;
     }
}

void G4ChipsKaonBuilder::attachProcess(G4ParticleDefinition * pDef)
{
     if ( verb > 0 ) {
        G4cout << " Using G4Qinelastic for " << pDef->GetParticleName()
	      << G4endl; 
     }
     G4ProcessManager* pMan = pDef->GetProcessManager();
     pMan->AddDiscreteProcess(theInelastic);
}
