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
// ClassName:   FTFP_BERT_TRV
//
// Author: 2009 John Apostolakis
//
//   created from FTFP_BERT - simple variant changing inelastic thresholds
//
// Modified:
// 16.08.2010 H.Kurashige: Remove inclusion of G4ParticleWithCuts 
// 18.10.2011 A.Ribon: Replace CHIPS with FTF/Preco for nuclear capture
//                     at rest of anti-protons.
// 27.07.2012 A.Ribon:  Use the new class G4BertiniAndFritiofStoppingPhysics
// 16.10.2012 A.Ribon:  Renamed the physics classes used
// 05.06.2014 A.Ribon:  Use (temporary) G4HadronHElasticPhysics
// 08.01.2015 A.Ribon:  Switch on low-mass diffraction dissociation
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "globals.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "G4HadronPhysicsFTFP_BERT_TRV.hh"
#include "FTFP_BERT_TRV.hh"

#include "G4WarnPLStatus.hh"

#include "G4EmParameters.hh"


FTFP_BERT_TRV::FTFP_BERT_TRV(G4int ver)
{
  // default cut value  (1.0mm) 
  // defaultCutValue = 1.0*CLHEP::mm;
  G4cout << "<<< Geant4 Physics List simulation engine: FTFP_BERT_TRV "<<G4endl;
  G4cout <<G4endl;
  defaultCutValue = 0.7*CLHEP::mm;  
  SetVerboseLevel(ver);

  G4WarnPLStatus exp;
  exp.Experimental("FTFP_BERT_TRV");

  // EM Physics
  G4EmStandardPhysicsGS* gsPhysics = new G4EmStandardPhysicsGS( ver );
  G4EmParameters::Instance()->SetMscStepLimitType( fUseSafety );
  RegisterPhysics( gsPhysics );

  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays 
  RegisterPhysics( new G4DecayPhysics(ver) );

   // Hadron Elastic scattering
  RegisterPhysics( new G4HadronHElasticPhysics(ver, true) );

   // Hadron Physics
  RegisterPhysics(  new G4HadronPhysicsFTFP_BERT_TRV(ver));

  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(ver) );

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver));
  
  // Neutron tracking cut
  RegisterPhysics( new G4NeutronTrackingCut(ver));

}

