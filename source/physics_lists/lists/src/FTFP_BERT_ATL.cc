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
// Author: Alberto Ribon
// Date:   April 2016
//
// New physics list FTFP_BERT_ATL.
// This is a modified version of the FTFP_BERT physics list for ATLAS.
// The physics list FTFP_BERT_ATL has the transition between Bertini (BERT)
// intra-nuclear cascade model and Fritiof (FTF) string model in the
// energy region [9, 12] GeV (instead of [3, 6] GeV as in FTFP_BERT).
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "globals.hh"
#include "G4ios.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "FTFP_BERT_ATL.hh"
#include "G4HadronPhysicsFTFP_BERT_ATL.hh"

#include "G4WarnPLStatus.hh"
#include "G4FTFTunings.hh"

FTFP_BERT_ATL::FTFP_BERT_ATL(G4int ver)
{
  if(ver > 0) {
    G4cout << "<<< Geant4 Physics List simulation engine: FTFP_BERT_ATL"<<G4endl;
    G4cout <<G4endl;
    G4WarnPLStatus exp;
    exp.Experimental("FTFP_BERT_ATL");
  }
  defaultCutValue = 0.7*CLHEP::mm;  
  SetVerboseLevel(ver);

  // Use the 4th tunes of Fritiof (FTF) string model, meant to to overcome
  // the problem of too optimistic (i.e. narrow) pion shower energy resolutions
  // in ATLAS calorimeters with respect to test-beam data.
  G4FTFTunings::Instance()->SetTuneApplicabilityState( 4, 1 );
  
  // EM Physics
  RegisterPhysics( new G4EmStandardPhysics(ver));

  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays 
  RegisterPhysics( new G4DecayPhysics(ver) );

   // Hadron Elastic scattering
  RegisterPhysics( new G4HadronElasticPhysics(ver) );

   // Hadron Physics
  RegisterPhysics( new G4HadronPhysicsFTFP_BERT_ATL(ver) );

  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(ver) );

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver));
  
  // Neutron tracking cut
  RegisterPhysics( new G4NeutronTrackingCut(ver));
}

