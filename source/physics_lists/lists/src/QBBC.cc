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
// $Id: QBBC.cc,v 1.4.2.1 2010/04/01 08:54:38 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-03-patch-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:QBBC
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
// 24.11.06 Add G4HadronHElasticPhysics and G4NeutronTrackingCut
// 16.05.2007 V.Ivanchenko: rename EM builders
//
//----------------------------------------------------------------------------
//

#include "QBBC.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmExtraPhysics.hh"
#include "G4QStoppingPhysics.hh"

#include "G4DataQuestionaire.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronDElasticPhysics.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4ChargeExchangePhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4NeutronTrackingCut.hh"

QBBC::QBBC( G4int ver, const G4String& type )
{
  G4DataQuestionaire it(photon, neutron);
  G4cout << "<<< Geant4 Physics List simulation engine: QBBC of type <"
	 << type <<">" <<G4endl;	
  G4cout <<G4endl;
  defaultCutValue = 0.7*mm;  
  SetVerboseLevel(ver);

  // EM Physics
  //RegisterPhysics( new G4EmStandardPhysics(ver));
  RegisterPhysics( new G4EmStandardPhysics_option2(ver));

  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics("extra EM"));

  // Decays
  RegisterPhysics( new G4DecayPhysics("decay",ver) );

   // Hadron Physics
  RegisterPhysics( new G4HadronHElasticPhysics(ver,false,type));
  //RegisterPhysics( new G4HadronElasticPhysics("elastic",ver,false));
  RegisterPhysics( new G4QStoppingPhysics("stopping",ver));
  RegisterPhysics( new G4IonBinaryCascadePhysics("ionBIC"));
  RegisterPhysics( new G4HadronInelasticQBBC(type, ver));

  // Neutron tracking cut
  RegisterPhysics( new G4NeutronTrackingCut("Neutron tracking cut", ver));
}		 

QBBC::~QBBC() 
{}

void QBBC::SetCuts()
{
  SetCutsWithDefault();   
}

